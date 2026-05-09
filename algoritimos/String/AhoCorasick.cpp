#include <bits/stdc++.h>
using namespace std;
//O(mk) m=tamanho total das strings, k=tamanho do alfabeto
// Chamar build() depois de adicionar todas as strings
const int K = 26;

struct Vertex {
    int next[K];
    int go[K];
    bool output = false;
    int link = 0;

    Vertex() {
        fill(begin(next), end(next), -1);
        fill(begin(go), end(go), 0);
    }
};

vector<Vertex> t(1);

void add_string(string const& s) {
    int v = 0;
    for (char ch : s) {
        int c = ch - 'a';
        if (t[v].next[c] == -1) {
            t[v].next[c] = t.size();
            t.emplace_back();
        }
        v = t[v].next[c];
    }
    t[v].output = true;
}

void build() {
    queue<int> q;
    for (int c = 0; c < K; c++) {
        if (t[0].next[c] == -1) {
            t[0].go[c] = 0;
        } else {
            t[0].go[c] = t[0].next[c];
            q.push(t[0].next[c]);
        }
    }
    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int c = 0; c < K; c++) {
            if (t[v].next[c] != -1) {
                int u = t[v].next[c];
                t[u].link = t[t[v].link].go[c];
                t[v].go[c] = u;
                q.push(u);
            } else {
                t[v].go[c] = t[t[v].link].go[c];
            }
        }
    }
}
