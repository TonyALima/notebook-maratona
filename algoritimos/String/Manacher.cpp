#include <bits/stdc++.h>
using namespace std;
//Funcao para Palindromos
const int S=200010;//Tamanho maximo da palavra

int odd[S];//Guarda o maior tamanho de palindromo centrado na posicao i, tamanho = 2*odd[i]-1
int even[S];//Guarda o maior tamanho de palindromo centrado nas posicoes i e i-1, tamanho = 2*even[i]

void manacher(string s) {//O(n)
    string t;
    for(auto c: s) {
        t += string("#") + c;
    }
    s = t + "#";
    int n = s.size();
    s = "$" + s + "^";
    vector<int> p(n + 2);
    int l = 0, r = 1;
    for(int i = 1; i <= n; i++) {
        p[i] = min(r - i, p[l + (r - i)]);
        while(s[i - p[i]] == s[i + p[i]]) {
            p[i]++;
        }
        if(i + p[i] > r) {
            l = i - p[i], r = i + p[i];
        }
    }
    for(int i=0;i<n;i++){
        if(i&1) odd[(i-1)/2]=p[i+1]/2;
        else even[i/2]=(p[i+1]-1)/2;
    }
}
