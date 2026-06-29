#include <bits/stdc++.h>
using namespace std;
// deteccao de ciclos, tempo de entrada e saida, nivel caso seja arvore

const int N=11;

int cnt=1;
int in[N],out[N],depth[N];//inicializa in e out com 0
int dis[N],low[N];
vector<int> adj[N];
set<int> artPoints;
set<pair<int,int>> bridges;


int DFS(int v,int nivel){//DFS com tempo de in e out, detecta ciclo e permite dizer se vertice faz parte da subarvore do outro
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

void dfs(int v,int par){//dfs com lowlink e dis para deteccao de pontes e pontos de articulacao
    // parentEdge ignora apenas UMA aresta de volta ao pai; em multigrafos arestas paralelas ao pai nao sao ignoradas
    bool parentEdge=false;
    int children=0;
    dis[v] = low[v] = cnt++;
    for(auto i:adj[v]){
        if(i==par&&!parentEdge){
            parentEdge=true;
            continue;
        }
        if(!dis[i]){
            children++;
            dfs(i,v);
            low[v] = min(low[v], low[i]);
            if (par!=-1 && low[i] >= dis[v]) artPoints.insert(v);
            if (low[i] > dis[v]) bridges.insert({min(v,i),max(v,i)});
        }else low[v] = min(low[v], dis[i]);
    }
    if (par==-1 && children>1) artPoints.insert(v);
}