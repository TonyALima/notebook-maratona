#include <bits/stdc++.h>
using namespace std;
//Usa BFS
//Necessario fonte e sumidouro

int fluxoMaximo(int ini, int fim, int tam){
    int u, v;
    int fluxo = 0;
    int bot;
    int rgrafo[tam][tam];
    for(u = 0; u < tam; u++)
        for(v = 0; v < tam; v++) rgrafo[u][v]=mat[u][v];
    
    while (BFS(ini, fim, tam)){
        bot = INF;
        for (v = fim; v != ini; v = anterior[v]){
            u = anterior[v];
            bot = min(bot, rgrafo[u][v]);
        }
        for (v = fim; v != ini; v = anterior[v]){
            u = anterior[v];
            rgrafo[u][v] -= bot; rgrafo[v][u] += bot; 
        }
        fluxo += bot;
    }
    return fluxo;
}
