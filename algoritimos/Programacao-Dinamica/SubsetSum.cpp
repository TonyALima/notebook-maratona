#include <bits/stdc++.h>
using namespace std;
//O(n*sum) nao funciona com numero negativo
vector<int> nums;

bool subsetsum(int sum) {
    vector<bool> prev(sum+1,false),cur(sum+1,false);
    prev[0]=true;
    for(auto i:nums){
        for(int j=0;j<=sum;j++){
            if(j<i) cur[j]=prev[j];
            else cur[j]=(prev[j]||prev[j-i]);
        }
        prev=cur;
    }
    return prev[sum];
}
