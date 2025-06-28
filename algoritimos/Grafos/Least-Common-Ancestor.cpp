#include <bits/stdc++.h>
using namespace std;

const int N=212345;
const int LN=20;
int depth[N],memo[LN][N];
vector<int> filhos[N];
//O(nlogn)
void dfs(int s){
    for(int i=0;i<filhos[s].size();i++){
        int f=filhos[s][i];
        if(memo[0][s]==f) continue;
        depth[f]=depth[s]+1;
        memo[0][f]=s;
        dfs(f);
    }
}
void compute(int n){
    for(int passos=1;passos<LN;passos++){
        for(int i=0;i<n;i++){
            memo[passos][i]=memo[passos-1][memo[passos-1][i]];
        }
    }
}
int walk(int s,int dis){
    int jump=0;
    while (dis){
        if(dis&1) s=memo[jump][s];
        jump++;
        dis>>=1;
    }
    return s;
}
int lca(int x,int y){
    if(depth[x]<depth[y]) swap(x,y);
    x=walk(x,depth[x]-depth[y]);
    if(x==y) return x;
    for(int k=LN-1;k>=0;k--){
        if(memo[k][x]!=memo[k][y]){
            x=memo[k][x];
            y=memo[k][y];
        }
    }
    return memo[0][x];
}
