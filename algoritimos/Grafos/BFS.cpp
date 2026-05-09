#include <bits/stdc++.h>
using namespace std;
//O(n+a) sendo a o numero de arestas
#define NMAX 1000
vector<int> adj[NMAX];
int vis[NMAX], anterior[NMAX];

int BFS(int ini, int fim, int tam){ // 1 se tiver caminho, 0 caso nao
    queue<int> fila;
    for (int i = 0; i < tam; i++) vis[i] = 0;
    fila.push(ini); vis[ini] = 1; anterior[ini] = -1;
    while (!fila.empty()){
        int u = fila.front(); fila.pop();
        for (int v : adj[u]){
            if (!vis[v]){
                if (v == fim){ anterior[v] = u; return 1; }
                fila.push(v); anterior[v] = u; vis[v] = 1;
            }
        }
    }
    return 0;
}
