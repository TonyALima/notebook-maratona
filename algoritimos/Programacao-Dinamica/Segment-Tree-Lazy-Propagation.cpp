#include <bits/stdc++.h>
using namespace std;
// range query/update, soma
const int maxn = 1000100;

int vet[maxn], tree[4 * maxn], lz[4 * maxn];

void build(int node, int l, int r)
{
    if (l == r)
    {
        tree[node] = vet[l];
        return;
    }
    int mid = (l + r) / 2;
    build(2 * node + 1, l, mid);
    build(2 * node + 2, mid + 1, r);
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

int main()
{
    int n = 5; // tamanho do vetor
    for (int i = 0; i < n; i++)
        vet[i] = 0;
    build(0, 0, n - 1);
    // Incrementa 3 no intervalo [1, 3]
    update(0, 0, n - 1, 1, 3, 3);
    // Soma total do intervalo [0, 4] -> deve ser 9 (0+3+3+3+0)
    cout << "Soma de [0, 4]: " << query(0, 0, n - 1, 0, 4) << '\n';
    return 0;
}
