#include <bits/stdc++.h>
using namespace std;

int mat[][];
int vis[], anterior[];

int BFS(int ini, int fim, int tam){ // 1 se tiver caminho, 0 caso nao
    int i;
    queue<int> fila;
    for (i = 0; i < tam; i++){ vis[i] = 0; }
    fila.push(ini); vis[ini] = 1; anterior[ini] = -1;
    while (!fila.empty()){
        for (i = 0; i < tam; i++){
            if (!vis[i] && mat[fila.front()][i]){
                if (i == fim) {anterior[i] = fila.front(); return 1;}
                fila.push(i); anterior[i] = fila.front(); vis[i] = 1;
            }
        }
        fila.pop();
    }
    return 0;
}
