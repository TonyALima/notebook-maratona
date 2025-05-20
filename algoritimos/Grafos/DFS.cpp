#include <bits/stdc++.h>
using namespace std;
// deteccao de ciclos

int mat[][];
int vis[], rec[]; // inicializa com 0

int DFS(int v, int tam){
    if(!vis[v]){
        vis[v] = 1;
        rec[v] = 1;
        for (int i = 0; i < tam; i++){
            if (mat[v][i]){
                if(!vis[i] && DFS(i, tam)) return 1;
                else if(rec[i]) return 1;
            }
        }
    }
    rec[v] = 0;
    return 0;
}
