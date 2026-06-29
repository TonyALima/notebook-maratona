#include <bits/stdc++.h>
using namespace std;
// range query/update, soma
//O(n)
const int maxn = 1000100;

struct SegTree {
    int tree[4 * maxn], lz[4 * maxn];

    void build(int* vet, int node, int l, int r)
    {
        lz[node] = 0;
        if (l == r)
        {
            tree[node] = vet[l];
            return;
        }
        int mid = (l + r) / 2;
        build(vet, 2 * node + 1, l, mid);
        build(vet, 2 * node + 2, mid + 1, r);
        tree[node] = tree[2 * node + 1] + tree[2 * node + 2];
    }

    void unlazy(int node, int tl, int tr)
    {
        if (lz[node] == 0)
            return;
        tree[node] += (tr - tl + 1) * lz[node];
        if (tl != tr)
        {
            lz[2 * node + 1] += lz[node];
            lz[2 * node + 2] += lz[node];
        }
        lz[node] = 0;
    }

    void update(int node, int tl, int tr, int l, int r, int v)
    {
        unlazy(node, tl, tr);
        if (tl > r || tr < l)
            return;
        if (tl >= l && tr <= r)
        {
            lz[node] += v;
            unlazy(node, tl, tr);
            return;
        }
        int mid = (tl + tr) / 2;
        update(2 * node + 1, tl, mid, l, r, v);
        update(2 * node + 2, mid + 1, tr, l, r, v);
        tree[node] = tree[2 * node + 1] + tree[2 * node + 2];
    }

    int query(int node, int tl, int tr, int l, int r)
    {
        unlazy(node, tl, tr);
        if (r < tl || l > tr)
            return 0;
        if (l <= tl && r >= tr)
            return tree[node];
        int mid = (tl + tr) / 2;
        return query(2 * node + 1, tl, mid, l, r) +
               query(2 * node + 2, mid + 1, tr, l, r);
    }
};

// Must be global — each instance is ~32 MB, stack cannot hold it
SegTree st1, st2;

int main()
{
    int n = 5; // tamanho do vetor
    int A[n], B[n];
    for (int i = 0; i < n; i++) A[i] = 0;
    for (int i = 0; i < n; i++) B[i] = 0;
    st1.build(A, 0, 0, n - 1);
    st2.build(B, 0, 0, n - 1);
    // Incrementa 3 no intervalo [1, 3] em st1
    st1.update(0, 0, n - 1, 1, 3, 3);
    // Soma total do intervalo [0, 4] de st1 -> deve ser 9 (0+3+3+3+0)
    cout << "st1 [0, 4]: " << st1.query(0, 0, n - 1, 0, 4) << '\n';
    // st2 nao foi alterada
    cout << "st2 [0, 4]: " << st2.query(0, 0, n - 1, 0, 4) << '\n';
    return 0;
}
