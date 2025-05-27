#include <bits/stdc++.h>
using namespace std;
//O(n)
long long int vet[412345];
int mem[412345];
void dp(int n){
    stack<int> s;
    for(int i=n-1;i>=0;i--){
        while(!s.empty()&&vet[i]<=vet[s.top()]) s.pop();
        if(s.empty()) mem[i]=-1;
        else mem[i]=s.top();
        s.push(i);
    }
}

int main(){
    int n;
    long long int k;
    cin>>n>>k;
    for(int i=0;i<n;i++){
        cin>>vet[i];
        vet[i+n]=vet[i]-k*(n+i);
        vet[i]-=k*i;
    }
    dp(2*n);
    cout<<mem[0]%n+1;
    for(int i=1;i<n;i++){
        cout<<" "<<mem[i]%n+1;
    }
    cout<<endl;
    return 0;
}
