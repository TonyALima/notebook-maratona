#include <bits/stdc++.h>
using namespace std;
#define INF 0x3F3F3F3F
#define NMAX 100
// Caminho minimo com aresta negativa, caminho percorrido

int m[NMAX][NMAX], custo[NMAX], anterior[NMAX];

// Retorna true se existir ciclo negativo
bool bellmanFord(int s, int n){
    int i, j, k;
    for (i = 0; i < n; i++){
        custo[i] = INF;
        anterior[i] = -1;
    }
    custo[s] = 0;

    for (k = 0; k < n; k++)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                if (m[i][j] != 0 && custo[i] != INF && custo[j] > custo[i] + m[i][j]) {
                    custo[j] = custo[i] + m[i][j];
                    anterior[j] = i;
                }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (m[i][j] != 0 && custo[i] != INF && custo[j] > custo[i] + m[i][j])
                return true; // Ciclo negativo
    
    return false;
}
