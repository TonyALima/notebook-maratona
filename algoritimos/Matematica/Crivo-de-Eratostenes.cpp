#include <bits/stdc++.h>
using namespace std;
// Cria um vetor de fatores primos / todos os primos ate N
//O(n)
const int N=11234567;
int fprimos[N]; // inicializa com 0
//p=primo !fprimos[p]||fprimos[p]==p
//p!=primo fprimos[p]&&fprimos[p]!=p
void crivo(){
    fprimos[0] = fprimos[1] = 1;
    for(int i = 2; i*i < N; i++){
        if(fprimos[i] == 0) fprimos[i] = i;
        for (int j = i*i; j < N; j+=i) fprimos[j] = fprimos[i];
    }
}

