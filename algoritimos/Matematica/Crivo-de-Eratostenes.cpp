#include <bits/stdc++.h>
using namespace std;
// Cria um vetor de fatores primos / todos os primos ate N
// fprimos[i] = menor fator primo de i (fprimos[primo] = primo, fprimos[0]=fprimos[1]=1)
//O(n)
#define MAXN 10000000
int fprimos[MAXN]; // inicializa com 0

void crivo(int n){
    fprimos[0] = fprimos[1] = 1;
    for(int i = 2; i < n; i++){
        if(fprimos[i] == 0){
            fprimos[i] = i;
        }
        if ((long long)i * i < n){
            for (int j = i*i; j < n; j+=i){
                if(fprimos[j] == 0) fprimos[j] = i;
            }
        }
    }
}
