#include <bits/stdc++.h>

using namespace std;

int LIS(vector<int>& arr)
{
    // Binary search approach
    int n = arr.size();
    vector<int> ans;
    ans.push_back(arr[0]);
    for(int i=1;i<arr.size();i++){
        if (arr[i] > ans.back()) {

            // If the current number is greater
            // than the last element of the answer
            // vector, it means we have found a
            // longer increasing subsequence.
            // Hence, we append the current number
            // to the answer vector.
            ans.push_back(arr[i]);
        }
        else {

            // If the current number is not
            // greater than the last element of
            // the answer vector, we perform
            // a binary search to find the smallest
            // element in the answer vector that
            // is greater than or equal to the
            // current number.

            // The lower_bound function returns
            // an iterator pointing to the first
            // element that is not less than
            // the current number.
            int low = lower_bound(ans.begin(), ans.end(),
                                  arr[i])
                      - ans.begin();

            // We update the element at the
            // found position with the current number.
            // By doing this, we are maintaining
            // a sorted order in the answer vector.
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