#include <bits/stdc++.h>
using namespace std;
#define MAXN 10000

int parent[MAXN];
int w[MAXN];

int find_set(int v)
{
    if (v == parent[v])
        return v;
    return parent[v] = find_set(parent[v]);
}

void make_set(int v)
{
    parent[v] = v;
    w[v] = 1;
}

void union_sets(int a, int b)
{
    a = find_set(a);
    b = find_set(b);
    if (a != b)
    {
        if (w[a] < w[b])
            swap(a, b);
        parent[b] = a;
        w[a] += w[b];
    }
}