#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

using namespace std;

using ll = long long;
using ull = unsigned long long;

ll INF = 0x3f3f3f3f3f3f3f3f;

#define endl '\n'
#define rep(i, a, b) for(int i = (a); i < (b); i++)
#define popcount(a) __builtin_popcount(a)
#define lzcnt(a) __builtin_ctl(a)
#define tzcnt(a) __builtin_ctz(a)

#ifndef ONLINE_JUDGE
    #define NDEBUG
#endif


// NUMBER THEORY ALGORITHMS

/**
 * @brief Finding GCD(a, b) in O(log min(a, b)); extended GCD
**/
ll GCD(ll a, ll b)
{
    while (b) {
        a %= b;
        swap(a, b);
    }
    return a;
}

ll GCD(ll a, ll b, ll &x, ll &y)
{
    if (a == 0) {
        x = 0; y = 1;
        return b;
    }
    ll x1, y1;
    ll g = GCD(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return g;
}

/**
 * @brief Modular multiplicative inverse using EEA in O(log N)
**/
ll ModMultInverse(ll k, ll mod)
{
    ll x, y;
    ll g = GCD(k, mod, x, y);
    assert(g != 1);
    return (x % mod + mod) % mod;
}

/**
 * @brief Checks if number is prime; does not handle negative and ones, O(sqrt(N))
**/
bool IsPrime(ll k)
{
    assert(k >= 2);
    for(ll i = 2, sq = 4; sq <= k; sq += i * 2 + 1, i++)
        if (k % i == 0)
            return 0;
    return 1;
}

/**
 * @brief Number of coprime numbers less than N
**/
ll Phi(ll k)
{
    ll res = k;
    assert(k >= 1);
    for(ll i = 2, sq = 4; sq <= k; sq += i * 2 + 1, i++)
        if (k % i == 0) {
            while (k % i == 0)
                k /= i;
            res -= res / i;
        }
    if (k > 1)
        res -= res / k;
    return res;
}


// DATA STRUCTURES



/**
 * @brief 0-indexed Fenwick tree
**/
template <typename T>
class FenwickTree {
    private:
        vector<T> __tree;
        T (&__add_func)(T t1, T t2);
        T __zerovalue;
        int __size;

    public:
        FenwickTree(int size, T zerovalue, T (&add_func)(T t1, T t2)) : __add_func(add_func)
        {
            __size = size;
            __tree.resize(size);
            fill(__tree.begin(), __tree.end(), zerovalue);
            __zerovalue = zerovalue;
        }

        void update(int pos, T val)
        {
            assert(pos >= 0);
            assert(pos < __size);
            pos++;
            for(; pos <= __size; pos += (pos & -pos))
                __tree[pos - 1] = __add_func(__tree[pos - 1], val);
            return;
        }

        T get(int pos)
        {
            assert(pos >= 0);
            assert(pos < __size);
            T ans = __zerovalue;
            pos++;
            for(; pos > 0; pos -= (pos & -pos))
                ans = __add_func(ans, __tree[pos - 1]);
            return ans;
        }
};


/**
 * @brief 0-indexed segment tree
**/
template<typename T>
class SegmentTree {
    private:
        int __size;
        vector<T> __tree;
        T (&__merge_function)(T t1, T t2);
        T __zerovalue;

        void __update(int l, int r, int t, int pos, T val)
        {
            if (l == r) {
                __tree[t] = val;
                return;
            }
            int lt = t * 2 + 1;
            int rt = t * 2 + 2;
            int m = (l + r) / 2;
            if (pos <= m)
                __update(l, m, lt, pos, val);
            else
                __update(m + 1, r, rt, pos, val);
            __tree[t] = __merge_function(__tree[lt], __tree[rt]);
            return;
        }

        T __get(int l, int r, int sl, int sr, int t)
        {
            if (l == sl && r == sr)
                return __tree[t];
            int lt = t * 2 + 1;
            int rt = t * 2 + 2;
            int m = (l + r) / 2;
            if (sr <= m)
                return __get(l, m, sl, sr, lt);
            else if (sl > m)
                return __get(m + 1, r, sl, sr, rt);
            else {
                T lres = __get(l, m, sl, m, lt);
                T rres = __get(m + 1, r, m + 1, sr, rt);
                return __merge_function(lres, rres);
            }
        }

    public:
        SegmentTree(int size, T zerovalue, T (&merge_function)(T t1, T t2)) : __merge_function(merge_function)
        {
            __size = size;
            __tree.resize(4 * size);
            fill(__tree.begin(), __tree.end(), zerovalue);
            __zerovalue = zerovalue;
        }

        void update(int pos, T val)
        {
            assert(pos >= 0);
            assert(pos < __size);
            __update(0, __size - 1, 0, pos, val);
            return;
        }
        
        T get(int l, int r)
        {
            assert(l >= 0);
            assert(r < __size);
            assert(l <= r);
            return __get(0, __size - 1, l, r, 0);
        }
};

/**
 * @brief DSU with amortization by default.
 * 
 * Numeration from 0.
**/
class DSU {
    private:
        vector<int> __lnk;
        vector<int> __sz;
        int __size;
        bool __amortized;

        int __find_nonamortized(int k)
        {
            if (__lnk[k] == k)
                return k;
            return find(__lnk[k]);
        }

        int __find_amortized(int k)
        {
            if (__lnk[k] == k)
                return k;
            return __lnk[k] = find(__lnk[k]);
        }

    public:
        void clear();

        DSU(int sz) {
            __lnk.resize(sz);
            __sz.resize(sz);
            __size = sz;
            clear();
            __amortized = 1;
        }

        DSU(int sz, bool amortized)
        {
            __lnk.resize(sz);
            __sz.resize(sz);
            __size = sz;
            clear();
            __amortized = amortized;
        }

        void clear()
        {
            for(int i = 0; i < __size; i++) {
                __lnk[i] = i;
                __sz[i] = 1;
            }
            return;
        }

        int find(int k)
        {
            assert(k >= 0);
            assert(k < __size);
            if (__amortized)
                return __find_amortized(k);
            else
                return __find_nonamortized(k);
        }

        bool same(int a, int b)
        {
            return find(a) == find(b);
        }

        void unite(int a, int b)
        {
            a = find(a);
            b = find(b);
            if (__sz[a] < __sz[b])
                swap(a, b);
            __lnk[b] = a;
            __sz[a] += __sz[b];
            return;
        }

        /**
         * Warning: works only if a -> b or b -> a
         **/
        void separate(int a, int b)
        {
            assert(!__amortized);
            if (find(a) == b)
                swap(a, b);
            __sz[a] -= __sz[b];
            __lnk[b] = b;
            return;
        }

        int size_of(int k)
        {
            assert(k >= 0);
            assert(k < __size);
            return __sz[k];
        }
};

const int N = 1e5;
const int MOD = 1e9+7;