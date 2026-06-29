#include <bits/stdc++.h>
using namespace std;
//Potenciação rapida
//O(logn) n sendo power
// Sem reducao modular; para potencia modular adicione % MOD nas duas multiplicacoes
long long fast_power(long long base, long long power) {
    long long result = 1;
    while(power > 0) {

        if(power % 2 == 1) { // Can also use (power & 1) to make code even faster
            result = (result*base) ;
        }
        base = (base * base);
        power = power / 2; // Can also use power >>= 1; to make code even faster
    }
    return result;
}
