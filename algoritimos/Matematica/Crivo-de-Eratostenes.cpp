#include <bits/stdc++.h>
using namespace std;
// Cria um vetor de fatores primos / todos os primos ate N
//O(n)
int fprimos[]; // inicializa com 0

void crivo(int n){
    fprimos[0] = fprimos[1] = 1;
    for(int i = 2; i < n; i++){
        if(fprimos[i] == 0){
            fprimos[i] = i;
        }
        if (i * i < n){
            for (int j = i*i; j < n; j+=i){
                fprimos[j] = i;
            }
        }
    }
}
