#include <bits/stdc++.h>
using namespace std;
// deteccao de ciclos, tempo de entrada e saida, nivel caso seja arvore

const int N=11;

int cnt=1;
int in[N],out[N],depth[N];//inicializa in e out com 0
vector<int> adj[N];


int DFS(int v,int nivel){
    in[v]=cnt++;
    depth[v]=nivel;
    for(int i=0;i<adj[v].size();i++){
        if(in[adj[v][i]]&&!out[adj[v][i]]) return 1;
        if(!in[adj[v][i]]){
            if(DFS(adj[v][i],nivel+1)) return 1;
        }
    }
    out[v]=cnt++;
    return 0;
}
