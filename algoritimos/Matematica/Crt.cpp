#include <bits/stdc++.h>
using namespace std;
struct Congruence {
    long long a, m;
};
//Calcular Chinese Remainder Theorem, usa inverso modular
long long mod_inv(long long a, long long m) {
    if (a <= 1) return a;
    return m - (mod_inv(m % a, a) * (m / a) % m);
}

long long chinese_remainder_theorem(vector<Congruence> const& congruences) {
    long long M = 1;
    for (auto const& congruence : congruences) {
        M *= congruence.m;
    }

    long long solution = 0;
    for (auto const& congruence : congruences) {
        long long a_i = congruence.a;
        long long M_i = M / congruence.m;
        long long N_i = mod_inv(M_i, congruence.m);
        solution = (solution + a_i * M_i % M * N_i) % M;
    }
    return solution;
}
