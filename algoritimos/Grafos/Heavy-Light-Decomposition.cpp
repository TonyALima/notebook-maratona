#include <bits/stdc++.h>
using namespace std;

const int MAXN = 100005;

// --- Segment Tree with Lazy Propagation (range assign, range sum) ---
int seg[4*MAXN], lz[4*MAXN];
bool dirty[4*MAXN];

void push_down(int node, int l, int r) {
    if (!dirty[node]) return;
    int mid = (l + r) / 2;
    seg[node] = (r - l + 1) * lz[node];//Range assign, += if range increment
    if (l != r) {
        for (int c : {2*node, 2*node+1}) { lz[c] = lz[node]; dirty[c] = true; }
    }
    dirty[node] = false;
}

void seg_update(int node, int l, int r, int ql, int qr, int val) {
    push_down(node, l, r);
    if (ql > r || qr < l) return;
    if (ql <= l && r <= qr) { lz[node] = val; dirty[node] = true; push_down(node, l, r); return; }
    int mid = (l + r) / 2;
    seg_update(2*node, l, mid, ql, qr, val);
    seg_update(2*node+1, mid+1, r, ql, qr, val);
    seg[node] = seg[2*node] + seg[2*node+1];//Sum
}

int seg_query(int node, int l, int r, int ql, int qr) {
    push_down(node, l, r);
    if (ql > r || qr < l) return 0;
    if (ql <= l && r <= qr) return seg[node];
    int mid = (l + r) / 2;
    return seg_query(2*node, l, mid, ql, qr) + seg_query(2*node+1, mid+1, r, ql, qr);//Sum
}

// --- HLD ---
int sz[MAXN], dep[MAXN], pai[MAXN], nxt[MAXN], tin[MAXN];
vector<int> g[MAXN];
int timer_hld = 0;

void dfs_sz(int v, int p) {
    sz[v] = 1; pai[v] = p;
    for (auto& u : g[v]) {
        if (u == p) continue;
        dep[u] = dep[v] + 1;
        dfs_sz(u, v);
        sz[v] += sz[u];
        if (sz[u] > sz[g[v][0]] || g[v][0] == p) swap(u, g[v][0]);
    }
}

void dfs_hld(int v, int p, int h) {
    nxt[v] = h; tin[v] = timer_hld++;
    for (auto u : g[v]) {
        if (u == p) continue;
        dfs_hld(u, v, u == g[v][0] ? h : u);
    }
}

// Path query (sum) from u to v
int query_path(int u, int v, int n) {
    int res = 0;
    while(nxt[u] != nxt[v]) {
        if (dep[nxt[u]] < dep[nxt[v]]) swap(u, v);
        res += seg_query(1, 0, n-1, tin[nxt[u]], tin[u]);
        u = pai[nxt[u]];
    }
    if (dep[u] > dep[v]) swap(u, v);
    res += seg_query(1, 0, n-1, tin[u], tin[v]);
    return res;
}

// Path update (range assign) from u to v
void update_path(int u, int v, int val, int n) {
    while(nxt[u] != nxt[v]) {
        if (dep[nxt[u]] < dep[nxt[v]]) swap(u, v);
        seg_update(1, 0, n-1, tin[nxt[u]], tin[u], val);
        u = pai[nxt[u]];
    }
    if (dep[u] > dep[v]) swap(u, v);
    seg_update(1, 0, n-1, tin[u], tin[v], val);
}

// Subtree query
int query_subtree(int v, int n) {
    return seg_query(1, 0, n-1, tin[v], tin[v] + sz[v] - 1);
}

//Range assign subtree
void update_subtree(int v,int val,int n){
    return seg_update(1,0,n-1,tin[v],tin[v]+sz[v]-1,val);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, q;
    cin >> n >> q;

    for (int i = 0; i < n-1; i++) {
        int u, v; cin >> u >> v;
        g[u].push_back(v);
        g[v].push_back(u);
    }

    dep[1] = 0;
    dfs_sz(1, 1);
    dfs_hld(1, 1, 1);

    while (q--) {
        int t, u, v; cin >> t >> u >> v;
        if (t == 1) update_path(u, v, /* val */ 1, n);
        else        cout << query_path(u, v, n) << '\n';
    }

    return 0;
}
