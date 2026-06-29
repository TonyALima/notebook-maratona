#include <bits/stdc++.h>
using namespace std;
//O(nlogn)
// dp[0] = -INF sentinela: evita acesso a dp[-1] quando a[v] e menor que todos os elementos
#define MAXN 100100
const int INF = 1e9;

int a[MAXN];
vector<int> tree[MAXN];
int r[MAXN];
int vis[MAXN];
vector<int> dp = []{vector<int> v(MAXN+1, INF); v[0]=INT_MIN; return v;}(); // dp[0]=sentinela; dp[1..] LIS
void dfs(int v)
{
    vis[v] = 1;
    int l = upper_bound(dp.begin()+1, dp.end(), a[v]) - dp.begin(); // 1-indexed: dp[1..] armazena a LIS
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
    r[v] = lower_bound(dp.begin()+1, dp.end(), INF) - dp.begin() - 1; // comprimento da LIS da raiz ao vértice v
    dp[l] = ant;
}
