#include <bits/stdc++.h>
using namespace std;
//Inverso Modular para divisão com MOD
//O(logMOD) precisa de gcd(n,MOD)=1 MOD precisa ser primo
const long long MOD=998244353;

long long mod_inv(long long num){return (num*fpow(num,MOD-2))%MOD;}

long long fpow(long long base, long long power) {
    long long result = 1;
    while(power > 0) {
        if(power&1) result = (result*base)%MOD ;
        base = (base * base)%MOD;
        power>>=1;
    }
    return result;
}
