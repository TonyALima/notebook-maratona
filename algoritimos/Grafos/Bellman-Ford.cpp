#include <bits/stdc++.h>
using namespace std;
#define INF 0x3F3F3F3F
#define NMAX 100
// Caminho minimo com aresta negativa
// O(V*E)
struct Edge { int u, v, w; };
vector<Edge> edges;
int custo[NMAX], anterior[NMAX];

// Retorna true se existir ciclo negativo
bool bellmanFord(int s, int n){
    for (int i = 0; i < n; i++){
        custo[i] = INF;
        anterior[i] = -1;
    }
    custo[s] = 0;

    for (int k = 0; k < n - 1; k++)
        for (auto& e : edges)
            if (custo[e.u] != INF && custo[e.v] > custo[e.u] + e.w){
                custo[e.v] = custo[e.u] + e.w;
                anterior[e.v] = e.u;
            }

    for (auto& e : edges)
        if (custo[e.u] != INF && custo[e.v] > custo[e.u] + e.w)
            return true; // Ciclo negativo

    return false;
}
