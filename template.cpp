#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")

using namespace std;

typedef long long ll;
typedef unsigned long long ull;
// typedef __int128_t int128;

const int N = 100;
long long  INF = 0x3f3f3f3f3f3f3f3f;
const int MOD = 1e9+7;

#define endl '\n'
#define rep(i, a, b) for(int i = (a); i < (b); i++)
#define popcount(a) __builtin_popcount(a)
#define lzcnt(a) __builtin_ctl(a)
#define tzcnt(a) __builtin_ctz(a)


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
    if (g != 1)
        throw invalid_argument("Inverse finding failed; check if base is prime");
    return (x % mod + mod) % mod;
}

/**
 * @brief Checks if number is prime; does not handle negative and ones, O(sqrt(N))
**/
bool IsPrime(ll k)
{
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

// GRAPH ALGORITHMS