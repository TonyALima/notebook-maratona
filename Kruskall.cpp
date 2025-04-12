#include <bits/stdc++.h>
using namespace std;
// arvore geradora minima
// usa union find

int mat[][];
int pai[], rnk[];

int find(int u){
    while (u != pai[u]) u = pai[u];
    return pai[u];
}

int merge(int u, int v){
    u = find(u); v = find(v);
    if(rnk[u] > rnk[v]) pai[v] = u; else pai[u] = v;
    if(rnk[u] == rnk[v]) rnk[v]++;
}

int kruskall(int tam){
    vector<tuple<int, int, int>> vec;
    vector<tuple<int, int, int>>::iterator it;
    int i, j, res = 0;
    for (i = 0; i < tam; i++){
        pai[i] = i; rnk[i] = 0;
        for (j = 0; j < tam; j++) vec.push_back(make_tuple(mat[i][j], i, j));
    }
    sort(vec.begin(), vec.end());
    for(it = vec.begin(); it != vec.end(); it++)
        if (find(get<1>(*it)) != find(get<2>(*it))){
            res += get<0>(*it); merge(get<1>(*it), get<1>(*it));
        }
    return res;
}
