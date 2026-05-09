#include <bits/stdc++.h>
using namespace std;
//O(n+m)
const int N = 112345;
vector<int> adj[N];
bool vis[N],in[N];
vector<int> ans;

bool dfs(int v) {
    in[v]=vis[v]=true;
    for (int u : adj[v]) {
        if (!vis[u]) {
            if(dfs(u)) return true;
        }else if(in[u]) return true;
    }
    in[v]=false;
    ans.push_back(v);
    return false;
}

bool topological_sort(int n) {//Ordena topologicamente e caso seja impossivel retorna false
    ans.clear();
    for (int i = 0; i < n; ++i) {
        if (!vis[i]) {
            if(dfs(i)) return false;
        }
    }
    reverse(ans.begin(), ans.end());
    return true;
}