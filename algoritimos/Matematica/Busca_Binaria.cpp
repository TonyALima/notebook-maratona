#include <bits/stdc++.h>
using namespace std;

const int INF=0x3f3f3f3f;

bool test(int n,int valor);//Funcao de teste no vetor

int bb(int n){//Busca Binaria que procura um valor que compra os requisitos da funcao teste
    int lmin=0,lmax=INF;
    int res=-1;
    while(lmin<=lmax){
        int mid=(lmin+lmax)/2;
        if(testa(n,mid)){
            res=mid;
            lmin=mid+1;
        }
        else lmax=mid-1;
    }
    return res;
}
