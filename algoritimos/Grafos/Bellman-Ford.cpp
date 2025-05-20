#include <bits/stdc++.h>
using namespace std;
// Caminho minimo com aresta negativa, caminho percorrido

int m[][], custo[], anterior[];

void bellmanFord(int s, int n){
    int i, j, k;
    for (i = 0; i < n; i++){
        custo[i] = INF;
        anterior[i] = -1;
    }
    custo[s] = 0;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                if (k != j && m[j][k] != 0 && custo[k] > custo[j]+m[j][k]){
                    custo[k] = custo[j]+m[j][k];
                    anterior[k] = j;
                }
}
