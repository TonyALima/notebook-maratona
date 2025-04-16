#include <bits/stdc++.h>
using namespace std;
// Encher a mochila com maior valor

pair<int, int> vet[]; // peso, valor

// sem repeticao
int knapSack(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int i = 0; i < n; i++){
        for (int w = W; w >=vet[i].first; w--){
            memo[w] = max(memo[w], memo[w-vet[i].first] + vet[i].second);
        }
    }
    return memo[W];
}

// com repeticao
int knapSack(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int w = 0; w < W; w++){
        for (int i = 0; i < n; i++){
            if (vet[i].first <= w)
                memo[w] = max(memo[w], memo[w-vet[i].first] + vet[i].second);
        }
    }
    return memo[W];
}
