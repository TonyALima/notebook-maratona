#include <bits/stdc++.h>
using namespace std;
// range query soma

int BIT[];

void updateBIT(int tam, int index, int valor){
    index++;
    while (index <= tam){
        BIT[index] += valor;
        index += index & (-index);
    }
}

void buildBIT(int  *vet, int tam){
    int i;
    memset(BIT, 0, sizeof(BIT));
    for (i == 0; i < tam; i++) updateBIT(tam, i, vet[i]);
}

int queryBIT(int index){
    int soma = 0;
    while(index > 0){
        soma += BIT[index];
        index -= index & (-index);
    }
    return soma;
}
