#include <bits/stdc++.h>
using namespace std;

// Requer: Template (Vetor, ccw)
// Se Vetor::x,y forem ll, troque "long double s=0" por "ll s=0" (soma exata) e
// "return fabsl(s)/2" por "return llabs(s)" retornando o dobro da area (evita .5 fracionario).
long double area(vector<Vetor>& poly){//Funciona para poligonos convexos ou nao
    long double s=0;
    int n=poly.size();
    for(int i=0;i<n;i++)
        s += poly[i] ^ poly[(i+1)%n];
    return fabsl(s)/2;
}

long double perimeter(vector<Vetor> &poly){
    long double s=0;
    int n=poly.size();
    for(int i=0;i<n;i++)
        s += ~(poly[i]-poly[(i+1)%n]);
    return s;
}

bool inPolygon(vector<Vetor> &polygon,Vetor p){//Poligono convexo em anti-horario, O(log n)
    Vetor piv=polygon[0];
    if(p==piv) return true;
    if(ccw(p,piv,polygon[1])<0) return false;
    if(ccw(polygon[polygon.size()-1],piv,p)<0) return false;
    int l=1,r=polygon.size()-1;
    while(l+1<r){
        int mid = (l+r)>>1;
        if(ccw(piv,polygon[mid],p)>=0) l=mid;
        else r=mid;
    }
    return ccw(polygon[l],polygon[r],p)>=0;
}
