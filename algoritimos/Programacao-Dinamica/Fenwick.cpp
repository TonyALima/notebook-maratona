#include <bits/stdc++.h>
using namespace std;
// range query soma
//O(logn)
#define NMAX 1000

struct Fenwick {
    int BIT[NMAX];

    void update(int tam, int index, int valor) {
        index++;
        while (index <= tam) {
            BIT[index] += valor;
            index += index & (-index);
        }
    }

    void build(int* vet, int tam) {
        memset(BIT, 0, sizeof(BIT));
        for (int i = 0; i < tam; i++) update(tam, i, vet[i]);
    }

    int query(int index) {
        int soma = 0;
        while (index > 0) {
            soma += BIT[index];
            index -= index & (-index);
        }
        return soma;
    }
};
