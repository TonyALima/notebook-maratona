#include <bits/stdc++.h>
using namespace std;
//Calcula a maior subsequencia crescente
//O(nlogn)
int LIS(vector<int>& arr)
{
    // Binary search approach
    int n = arr.size();
    vector<int> ans;
    ans.push_back(arr[0]);
    for(int i=1;i<arr.size();i++){
        if (arr[i] > ans.back()) {
            // Se numero foi maior bota ele na sequencia
            ans.push_back(arr[i]);
        }
        else {
            //Senao troca o menor numero maior q ele por ele
            int low = lower_bound(ans.begin(), ans.end(),
                                  arr[i])
                      - ans.begin();
            ans[low] = arr[i];
        }
    }
    return ans.size();
}

int main(){
    int n;
    vector<int> v;
    cin>>n;
    while(n--){
        int aux;
        cin>>aux;
        v.push_back(aux);
    }
    cout<<LIS(v)<<endl;
    return 0;
}
