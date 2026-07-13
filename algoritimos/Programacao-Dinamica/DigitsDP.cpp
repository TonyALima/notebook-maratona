#include <bits/stdc++.h>
using namespace std;

// Generic digit DP template. Edit the 3 pieces marked EDIT to change what
// "characteristic" of a number is being counted.

const int N = 80000; // max digits; recursion depth (not compute time) is the real limit here

// ---- 1) EDIT: State (here: digit sum mod 9) ----
struct State {
    int sum;
};
const int STATE_SIZE = 9;//Depends on the State
int encode(const State& s) { return s.sum; }
State initialState() { return State{0}; }

// ---- 2) EDIT: transition (how digit d updates the state) ----
State transition(const State& s, int d, bool started) {
    State ns = s;
    ns.sum = (s.sum + d) % 9;
    return ns;
}

// ---- 3) EDIT: isValid (which final states count; started = leading zeros done) ----
long long isValid(const State& s, bool started) {
    if (!started) return 1; // the number 0 itself
    return s.sum == 0;      // example: digit sum divisible by 9
}

// ---- DP machinery, rarely needs changes ----

vector<int> digits;
long long memo[N][2][STATE_SIZE];
bool visited[N][2][STATE_SIZE];

long long solve(int pos, bool tight, bool started, State state) {
    if (pos == (int)digits.size()) {
        return isValid(state, started);
    }

    int stCode = encode(state);
    if (!tight && visited[pos][started][stCode]) {
        return memo[pos][started][stCode];
    }

    int limit = tight ? digits[pos] : 9;
    long long res = 0;
    for (int d = 0; d <= limit; d++) {
        bool nStarted = started || d != 0;
        State nState = transition(state, d, nStarted);
        res += solve(pos + 1, tight && d == limit, nStarted, nState);
    }

    if (!tight) {
        visited[pos][started][stCode] = true;
        memo[pos][started][stCode] = res;
    }
    return res;
}

long long countUpTo(long long n) {
    if (n < 0) return 0;
    digits.clear();
    if (n == 0) digits.push_back(0);
    for (long long x = n; x > 0; x /= 10) digits.push_back(x % 10);
    reverse(digits.begin(), digits.end());
    memset(visited, 0, sizeof(visited));
    return solve(0, true, false, initialState());
}

long long solve_range(long long l, long long r) {
    return countUpTo(r) - countUpTo(l - 1);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    long long l, r;
    cin >> l >> r;
    cout << solve_range(l, r) << "\n";
    return 0;
}
