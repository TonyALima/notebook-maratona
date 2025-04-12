#include <bits/stdc++.h>
using namespace std;

typedef struct{double x, y;}Ponto;

// Função para calcular o produto vetorial
double crossProduct(Ponto a, Ponto b) {
    return a.x * b.y - a.y * b.x;
}

bool doIntersect(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    double d1 = crossProduct({p3.x - p1.x, p3.y - p1.y}, {p2.x - p1.x, p2.y - p1.y});
    double d2 = crossProduct({p4.x - p1.x, p4.y - p1.y}, {p2.x - p1.x, p2.y - p1.y});
    double d3 = crossProduct({p1.x - p3.x, p1.y - p3.y}, {p4.x - p3.x, p4.y - p3.y});
    double d4 = crossProduct({p2.x - p3.x, p2.y - p3.y}, {p4.x - p3.x, p4.y - p3.y});

    return (d1 * d2 < 0) && (d3 * d4 < 0);
}
