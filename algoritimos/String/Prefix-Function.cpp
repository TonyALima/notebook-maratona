#include <bits/stdc++.h>
using namespace std;

const int N=1123456;
int pi[N];//Caso for acessar tomar cuidado pois o valor dentro eh o tamanho
int ans[N];//Prefixo de tamanho i aparece ans[i] vezes
//Usado para contar prefixo, procurar substring e contar substrings
void prefix_function(string s){//Guarda o tamanho da string que eh sufixo e prefixo simultaneamente
    //A letra i eh o final de um sufixo e prefixo de tamanho pi[i]
    pi[0]=0;
    for(int i=1;i<s.length();i++){
        int j=pi[i-1];
        while(j>0&&s[j]!=s[i]) j=pi[j-1];
        if(s[j]==s[i]) pi[i]=j+1;
        else pi[i]=0; 
    }
}

void conta_prefixos(string s){//Conta quantas vezes cada prefixo aparece na string
    int n=s.length();
    prefix_function(s);
    for (int i = 0; i < n; i++)
        ans[pi[i]]++;
    for (int i = n-1; i > 0; i--)
        ans[pi[i-1]] += ans[i];
    for (int i = 0; i <= n; i++)
        ans[i]++;
}
