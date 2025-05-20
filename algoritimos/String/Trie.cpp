#include <bits/stdc++.h>
using namespace std;
//Cria um dicionario de prefixos
const int K = 26;
struct Vertex {
    int next[K];
    int output = 0;
    bool eliminado = false;
    Vertex() {
        fill(begin(next), end(next), -1);
    }
};
vector<Vertex> trie(1);
//Adiciona o prefixo e da update no output
int add_string(string const& s) {
    int v = 0;
    for (char ch : s) {
        int c = ch - 'a';
        if(trie[v].eliminado) return 0;
        if (trie[v].next[c] == -1) {
            trie[v].next[c] = trie.size();
            trie.emplace_back();
        }
        v = trie[v].next[c];
    }
    if(trie[v].eliminado) return 0;
    trie[v].output++;
    return 1;
}
//Conta quantos prefixos terminam dps do ponto
int conta_out(int v){
    int res=0;
    if(trie[v].eliminado) return 0;
    for(int i=0;i<K;i++){
        if(trie[v].next[i]!=-1){
            res+=conta_out(trie[v].next[i]);
        }
    }
    res+=trie[v].output;
    return res;
}
//Apaga o prefixo e conta os eliminados
int remove_string(string const& s){
    int v = 0;
    for (char ch : s) {
        int c = ch - 'a';
        if(trie[v].eliminado) return 0;
        if (trie[v].next[c] == -1) {
            trie[v].next[c] = trie.size();
            trie.emplace_back();
        }
        v = trie[v].next[c];
    }
    int res=conta_out(v);
    trie[v].eliminado = true;
    return res;
}


int main(){
    int n,res=0;
    cin>>n;
    while(n--){
        int tipo;
        string s;
        cin>>tipo;
        cin>>s;
        if(tipo==2) res+=add_string(s);
        else res-=remove_string(s);
        cout<<res<<endl;
    }
    return 0;
}
