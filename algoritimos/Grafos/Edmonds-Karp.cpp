#include <bits/stdc++.h>
using namespace std;
//O(Ef) E=num de arestas, f=fluxomax
struct Edge{
    int to;
    int cap;
};

const int N=112,INF=0x3f3f3f3f;

int vis[N],anterior[N],aresta[N];
vector<int> g[N];
vector<Edge> edges;
set<int> s;
//Funcao auxiliar para criar arestas
void addEdge(int a,int b,int c){
    g[a].push_back(edges.size());
    edges.push_back(Edge{b,c});
    g[b].push_back(edges.size());
    edges.push_back(Edge{a,0}); // aresta reversa com cap 0 (grafo direcionado); para nao-direcionado chame addEdge(b,a,c) separadamente
}
//Roda BFS para pegar o caminho
int BFS(int ini, int fim){ // 1 se tiver caminho, 0 caso nao
    int i;
    queue<int> fila;
    memset(vis,0,sizeof(vis));
    fila.push(ini); vis[ini] = 1; anterior[ini] = -1;
    while (!fila.empty()){
        int a=fila.front();fila.pop();
        for(auto i:g[a]){
            int b=edges[i].to,c=edges[i].cap;
            if(!vis[b]&&c){
                anterior[b]=a;
                aresta[b]=i;
                if(b==fim) return 1;
                fila.push(b);
                vis[b]=1;
            }
        }
    }
    return 0;
}

int fluxoMaximo(int ini, int fim){
    int u, v;
    int fluxo = 0;
    int bot;
    while (BFS(ini, fim)){
        bot = INF;
        for (v = fim; v != ini; v = anterior[v]){
            u = aresta[v];
            bot = min(bot, edges[u].cap);
        }
        for (v = fim; v != ini; v = anterior[v]){
            u = aresta[v];
            edges[u].cap -= bot; edges[u^1].cap += bot; 
        }
        fluxo += bot;
    }
    return fluxo;
}
//Funcao que determina arestas que participam do minCut
int minCut(int tam){
    for(int i=0;i<tam;i++){
        if(vis[i]){
            for(auto j:g[i]){
                int b=edges[j].to;
                if(!vis[b]){
                    s.insert(j/2+1);
                }
            }
        }
    }
    return s.size();
}
