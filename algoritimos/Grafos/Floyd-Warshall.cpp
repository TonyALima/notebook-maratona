#include <bits/stdc++.h>
using namespace std;
// Deteccao de ciclo negativo e caminho minimo para qualquer u, v

int m[][], custo[][];

bool floydWarshall(int n){
    int i, j, k;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i == j) custo[i][j] = 0;
            else custo[i][j] = m[i][j];
        }
    }

    for (k = 0; k < n; k++)
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
                if (custo[i][k] != INF && custo[k][j] != INF &&
                    custo[i][j] > custo[i][k] + custo[k][j])
                    custo[i][j] = custo[i][k] + custo[k][j];
            }
            if (custo[j][j]< 0)
                return true;
        }
    return false;
}
