#include <bits/stdc++.h>
using namespace std;

// Kuhn's Algorithm — Maximum Bipartite Matching
// Problem: given two groups L (left) and R (right) with edges between them,
// find the largest set of edges where every vertex appears at most once.

const int NONE = -1;

int L, R;                       // number of left and right vertices
vector<int> adj[505];           // adj[u] = list of right vertices u can match to
int matchL[505], matchR[505];   // matchL[u] = right vertex matched to u, and vice versa

// DFS from left vertex u looking for an augmenting path.
// Returns true if an augmenting path was found (and the matching was updated).
bool tryAugment(int u, vector<bool>& visited) {
    for (int v : adj[u]) {
        if (visited[v]) continue;   // already tried to re-route through v in this DFS
        visited[v] = true;

        // v is unmatched — we can directly match u to v
        // v is matched to matchR[v] — ask matchR[v] if it can go elsewhere
        if (matchR[v] == NONE || tryAugment(matchR[v], visited)) {
            // augmenting path found: update the matching
            matchL[u] = v;
            matchR[v] = u;
            return true;
        }
    }
    return false;   // no augmenting path from u
}

int maxMatching() {
    fill(matchL, matchL + L, NONE);
    fill(matchR, matchR + R, NONE);

    int result = 0;
    for (int u = 0; u < L; u++) {
        // visited is reset for each left vertex so that each DFS call
        // gets a fresh chance to reroute right vertices
        vector<bool> visited(R, false);
        if (tryAugment(u, visited))
            result++;
    }
    return result;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    L = 3; R = 3;
    adj[0] = {0, 1};
    adj[1] = {1};
    adj[2] = {1, 2};

    cout << "Maximum matching: " << maxMatching() << "\n";
    for (int u = 0; u < L; u++)
        if (matchL[u] != NONE)
            cout << "  Worker " << u << " -> Job " << matchL[u] << "\n";

    return 0;
}
