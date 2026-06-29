#include <bits/stdc++.h>
using namespace std;

// Requer: Template (Vetor, EPS, ccw)
// Se Vetor::x,y forem ll, as margens +-EPS abaixo sao desnecessarias (comparacoes ficam exatas).
bool onSegment(Vetor a, Vetor b, Vetor p) {//p pertence ao segmento ab
    return ccw(a,b,p)==0 &&
       min(a.x,b.x)-EPS <= p.x && p.x <= max(a.x,b.x)+EPS &&
       min(a.y,b.y)-EPS <= p.y && p.y <= max(a.y,b.y)+EPS;
}

bool intersect(Vetor a, Vetor b, Vetor c, Vetor d) {//segmentos ab e cd se intersectam
    int d1 = ccw(a,b,c), d2 = ccw(a,b,d);
    int d3 = ccw(c,d,a), d4 = ccw(c,d,b);
    if(d1*d2<0 && d3*d4<0) return true;
    if(d1==0 && onSegment(a,b,c)) return true;
    if(d2==0 && onSegment(a,b,d)) return true;
    if(d3==0 && onSegment(c,d,a)) return true;
    if(d4==0 && onSegment(c,d,b)) return true;
    return false;
}

//Segmento-Segmento
long double distSS(Vetor a,Vetor b,Vetor c,Vetor d){
    if(intersect(a,b,c,d)) return 0;
    return min({a.dist(c,d),b.dist(c,d),c.dist(a,b),d.dist(a,b)});
}
