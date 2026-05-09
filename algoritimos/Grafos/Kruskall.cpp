#include <bits/stdc++.h>
using namespace std;
// arvore geradora minima
// usa union find
//O(elogv)
#define NMAX 1000
int pai[NMAX], rnk[NMAX];

int find(int u){
    return pai[u] = (pai[u] == u ? u : find(pai[u]));
}

void merge(int u, int v){
    u = find(u); v = find(v);
    if(rnk[u] > rnk[v]) pai[v] = u; else pai[u] = v;
    if(rnk[u] == rnk[v]) rnk[v]++;
}

// edges: vetor de (peso, u, v)
int kruskall(int n, vector<tuple<int,int,int>>& edges){
    for (int i = 0; i < n; i++){
        pai[i] = i; rnk[i] = 0;
    }
    sort(edges.begin(), edges.end());
    int res = 0;
    for (auto& [w, u, v] : edges)
        if (find(u) != find(v)){
            res += w; merge(u, v);
        }
    return res;
}
