#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-6;
// Busca ternária para encontrar o ponto de máximo (ou mínimo, dependendo da função)
//O(logn) n = tamanho de busca
double ternary_search(std::function<double(double)> f, double left, double right) {
    while (right - left > EPS) {
        double m1 = left + (right - left) / 3.0;
        double m2 = right - (right - left) / 3.0;
        if (f(m1) < f(m2))  // Para máximo
            left = m1;
        else
            right = m2;
    }
    return (left + right) / 2.0;
}

/*Para numeros inteiros
int ternary_search(std::function<int(int)> f, int left, int right) {
    while (left<right) {
        int m1 = left + (right - left) / 3;
        int m2 = right - (right - left) / 3;
        if (f(m1) > f(m2))  // Para min
            left = m1+1;
        else
            right = m2-1;
    }
    return left;
}
*/

int main() {
    // Exemplo de função unimodal: f(x) = - (x - 2)^2 + 4, máximo em x = 2
    auto f = [](double x) {
        return -(x - 2) * (x - 2) + 4;
    };

    double a = 0.0, b = 4.0;//Limites superior e inferior

    // Encontra o ponto de máximo ou
    double xm = ternary_search(f, a, b);

    return 0;
}
