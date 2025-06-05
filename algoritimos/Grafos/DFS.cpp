#include <bits/stdc++.h>
using namespace std;
// deteccao de ciclos

int cnt=0;
int in[N],out[N],depth[N];//inicializa in e out com -1
vector<int> adj[N];


int DFS(int v,int nivel){
    int res=0;
    in[v]=cnt++;
    depth[v]=nivel;
    for(int i=0;i<adj[v].size();i++){
        if(in[adj[v][i]]!=-1&&out[adj[v][i]]==-1) return 1;
        if(in[adj[v][i]]==-1){
            if(DFS(adj[v][i],nivel+1)) res=1;
        }
    }
    out[v]=cnt++;
    return res;
}
