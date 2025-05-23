#include <bits/stdc++.h>
using namespace std;

#define MAXN 100100
const int INF = 1e9;

int a[MAXN];
vector<int> tree[MAXN];
int r[MAXN];
int vis[MAXN];
vector<int> dp(MAXN, INF);
void dfs(int v)
{
    vis[v] = 1;
    int l = upper_bound(dp.begin(), dp.end(), a[v]) - dp.begin(); // índice do primeiro elemento em dp maior que a[i]
    int ant = dp[l];
    if (dp[l - 1] < a[v] && a[v] < dp[l])
    {
        dp[l] = a[v];
    }
    for (auto u : tree[v])
    {
        if (!vis[u])
        {
            dfs(u);
        }
    }
    r[v] = lower_bound(dp.begin(), dp.end(), INF) - dp.begin() - 1; // comprimento da LIS da raiz ao vértice v
    dp[l] = ant;
}