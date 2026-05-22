#include <bits/stdc++.h>
using namespace std;

// Requer: Template (Vetor, ccw)
double area(vector<Vetor>& poly){//Funciona para poligonos convexos ou nao
    double s=0;
    int n=poly.size();
    for(int i=0;i<n;i++)
        s += poly[i] ^ poly[(i+1)%n];
    return fabs(s)/2;
}

bool inPolygon(vector<Vetor> polygon,Vetor p){//Poligono convexo em anti-horario, O(log n)
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
