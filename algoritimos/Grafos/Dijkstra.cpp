#include <bits/stdc++.h>
using namespace std;
// Menor caminho
//O(e+nlogn) sendo e=arestas
#define SIZE 1000
#define INF 0x3f3f3f3f
typedef pair<int, int> ii;
typedef vector<ii> vii;
typedef vector<int> vi;

vii adj[SIZE];
vi custo(SIZE, INF); // resetar com fill(custo.begin(), custo.end(), INF) entre chamadas

void dijkstra(int s){
    custo[s] = 0;
    priority_queue<ii, vii, greater<ii>> heap;
    heap.push(ii(0, s));
    while(!heap.empty()){
        ii menor = heap.top(); heap.pop();
        int d = menor.first, u = menor.second;
        if (d > custo[u]) continue;
        for (int i = 0; i < adj[u].size(); i++){
            ii v = adj[u][i];
            if(custo[u]+v.second < custo[v.first]){
                custo[v.first] = custo[u] + v.second;
                heap.push(ii(custo[v.first], v.first));
            }
        }
    }
}
