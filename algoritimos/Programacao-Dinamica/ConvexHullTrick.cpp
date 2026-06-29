#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

// O(n) or O(n log n) DP optimization
// Optimizes DP transitions of the form:
//   dp[i] = min over j of:  dp[j] + a[j] * x[i] + b[j]
// Each j defines a line:  y = a[j] * X + b[j]
// For a fixed query X = x[i], we want the minimum y among all lines.
// Only lines on the "lower envelope" can ever be optimal.
// The lower envelope, read left to right, is ordered largest slope → smallest slope.
// When slopes are added in DECREASING order AND queries are in INCREASING order,
// the optimal line index only moves forward → amortized O(1) per operation.
// For arbitrary query order, use the binary-search variant: O(log n) per query.

struct Line {
    ll m, b;
    ll eval(ll x) const { return m * x + b; }
};

// Returns true if line B is made redundant by lines A (left) and C (right).
bool bad(Line A, Line B, Line C) {
    return (__int128)(C.b - A.b) * (A.m - B.m) <= (__int128)(B.b - A.b) * (A.m - C.m);
}


// REQUIREMENT: addLine must be called with slopes in DECREASING order.
// REQUIREMENT: query must be called with x values in NON-DECREASING order.
//   hull[0] has the steepest slope → optimal for the smallest x values.
//   As x grows, the optimal line shifts rightward in the hull.
struct CHT_Monotone {
    vector<Line> hull;
    int ptr = 0;  // points to the current best line; only moves right

    void addLine(ll slope, ll intercept) {
        Line L = {slope, intercept};
        // Remove the last hull line if C makes it redundant.
        while (hull.size() >= 2 && bad(hull[hull.size()-2], hull[hull.size()-1], L))
            hull.pop_back();
        hull.push_back(L);
    }

    // x must be non-decreasing across successive calls.
    ll query(ll x) {
        // Once hull[ptr] is no longer beaten by hull[ptr+1], we stop
        while (ptr + 1 < (int)hull.size() && hull[ptr].eval(x) >= hull[ptr+1].eval(x))
            ptr++;
        return hull[ptr].eval(x);
    }
};

// Same envelope construction, but instead of a moving pointer we binary search.
// Slopes must still be added in DECREASING order; queries can be in any order.
struct CHT_BinarySearch {
    vector<Line> hull;

    void addLine(ll slope, ll intercept) {
        Line L = {slope, intercept};
        while (hull.size() >= 2 && bad(hull[hull.size()-2], hull[hull.size()-1], L))
            hull.pop_back();
        hull.push_back(L);
    }

    ll query(ll x) {
        int lo = 0, hi = (int)hull.size() - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            // If the next line is already better at x, the optimal is to the right
            if (hull[mid].eval(x) >= hull[mid+1].eval(x))
                lo = mid + 1;
            else
                hi = mid;
        }
        return hull[lo].eval(x);
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    // --- Simple demo: minimum of a set of lines ---

    CHT_Monotone cht;
    cht.addLine(3, 1);    // slope 3  (steepest, good for small x)
    cht.addLine(1, 5);    // slope 1
    cht.addLine(-1, 9);   // slope -1 (shallowest, good for large x)

    cout << "CHT Monotone (min of 3 lines, queries in increasing x):\n";
    cout << "x=0: " << cht.query(0) << " (expected 1)\n";
    cout << "x=2: " << cht.query(2) << " (expected 7)\n";
    cout << "x=4: " << cht.query(4) << " (expected 5)\n";
    cout << "x=6: " << cht.query(6) << " (expected 3)\n";

    // Same demo with binary search (queries can come in any order)
    CHT_BinarySearch cht2;
    cht2.addLine(3, 1);
    cht2.addLine(1, 5);
    cht2.addLine(-1, 9);

    cout << "\nCHT Binary Search (same lines, arbitrary query order):\n";
    cout << "x=6: " << cht2.query(6) << " (expected 3)\n";
    cout << "x=0: " << cht2.query(0) << " (expected 1)\n";
    cout << "x=4: " << cht2.query(4) << " (expected 5)\n";
    return 0;
}
