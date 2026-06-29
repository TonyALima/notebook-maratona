#include <bits/stdc++.h>
using namespace std;
//O(nlogn) construcao O(1)ouO(logn) para querys
const int K=21;
const int N=112345;

// Must be global — each instance is ~10 MB, stack cannot hold it
struct SparseTable {
    int st[K + 1][N];
    int lg[N + 1];

    int f(int a, int b) { return min(a, b); }

    void build(int* v, int n) {
        copy(v, v + n, st[0]);
        lg[1] = 0;
        for (int i = 2; i <= n; i++)
            lg[i] = lg[i/2] + 1;
        for (int i = 1; i <= K; i++)
            for (int j = 0; j + (1 << i) <= n; j++)
                st[i][j] = f(st[i-1][j], st[i-1][j + (1 << (i-1))]);
    }

    int queryI(int L, int R) { // Funcao Idempotente: O(1)
        int i = lg[R - L + 1];
        return f(st[i][L], st[i][R - (1 << i) + 1]);
    }

    long long query(int L, int R) { // Exemplo para soma: O(logn)
        long long sum = 0;
        for (int i = K; i >= 0; i--) {
            if ((1 << i) <= R - L + 1) {
                sum += st[i][L];
                L += 1 << i;
            }
        }
        return sum;
    }
};
