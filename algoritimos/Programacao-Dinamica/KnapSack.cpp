#include <bits/stdc++.h>
using namespace std;
// Encher a mochila com maior valor
//O(n*W)
int peso[], valor[];

// com repeticao
int knapSack(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int w = 0; w < W; w++){
        for (int i = 0; i < n; i++){
            if (peso[i] <= w)
                memo[w] = max(memo[w], memo[w-peso[i]] + valor[i]);
        }
    }
    return memo[W];
}

// sem repeticao
int knapSack(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int i = 0; i < n; i++){
        for (int w = W; w >=peso[i]; w--){
            memo[w] = max(memo[w], memo[w-peso[i]] + valor[i]);
        }
    }
    return memo[W];
}

// versao com matriz.

int mat[][];

int knapSack(int W, int n) {
    for (int i = 1; i <= n; ++i) {
        for (int w = 0; w <= W; ++w) {
            if (peso[i-1] <= w) {
                mat[i][w] = max(mat[i-1][w], mat[i-1][w - peso[i-1]] + valor[i-1]);
            } else {
                mat[i][w] = mat[i-1][w];
            }
        }
    }
    return mat[n][W];
}

vector<int> escolhidos(int W, int n){
    int w = W;
    vector<int> itens = vector<int>();
    while(n > 0 && w > 0){
        if (mat[n][w] != mat[n-1][w]) {
            itens.push_back(n-1);
            w -= peso[n-1];
        }
        n--;
    }
    return itens;
}
