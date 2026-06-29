#include <bits/stdc++.h>
using namespace std;

// Sqrt Decomposition for range queries with point updates.
// Core idea: divide the array into blocks of size ~sqrt(n).
// Query [l, r]:
//   - Left partial block  : iterate element by element
//   - Full blocks in middle: use precomputed block sum directly
//   - Right partial block  : iterate element by element
//   Cost: O(sqrt(n))
// Update index i:
//   - Update the element and adjust its block sum by the delta
//   Cost: O(1)
// Better than segment tree when the operation is hard to merge (distinct count, median)

struct SqrtDecomp {
    int n, B;           // array size, block size
    vector<long long> a;       // original array
    vector<long long> block;   // block[i] = sum of elements in block i

    SqrtDecomp(vector<int>& arr) {
        n = arr.size();
        B = max(1, (int)sqrt(n));   // block size ~ sqrt(n)
        a.assign(arr.begin(), arr.end());
        block.assign((n + B - 1) / B, 0);

        // Precompute block sums in O(n)
        for (int i = 0; i < n; i++)
            block[i / B] += a[i];
    }

    // Update a[i] to val in O(1)
    void update(int i, long long val) {
        block[i / B] += val - a[i];  // adjust block sum by the delta
        a[i] = val;
    }

    // Sum of a[l..r] in O(sqrt(n))
    long long query(int l, int r) {
        long long sum = 0;
        int lb = l / B;   // block containing l
        int rb = r / B;   // block containing r

        if (lb == rb) {
            // l and r are in the same block — just iterate
            for (int i = l; i <= r; i++)
                sum += a[i];
            return sum;
        }

        // Left partial block: from l to the end of block lb
        for (int i = l; i < (lb + 1) * B; i++)
            sum += a[i];

        // Full blocks in between: O(n/B) = O(sqrt(n)) iterations
        for (int b = lb + 1; b < rb; b++)
            sum += block[b];

        // Right partial block: from start of block rb to r
        for (int i = rb * B; i <= r; i++)
            sum += a[i];

        return sum;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    vector<int> arr = {1, 3, 2, 7, 4, 6, 5, 8, 2, 3};
    SqrtDecomp sd(arr);
    sd.query(2, 7);//Query from position 2 to 7
    sd.update(3, 10);//Update element in position 3 to 10
    sd.query(2, 7);
    sd.query(0, 9);

    return 0;
}
