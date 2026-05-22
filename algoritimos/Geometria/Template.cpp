#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const double EPS = 1e-9;
const double PI = acos(-1);

// Vetor: +(add) -(sub) *(double)=escala *(Vetor)=dot ^=cross ~=modulo dist(p,q)=dist ao segmento pq
// cmpPolar: usar com sort + Vetor::pivot para ordenar por angulo polar
class Vetor{
    public:
    double x,y;//Pode ser trocado por long long int dependendo da questao
    Vetor operator+(Vetor q) const{
        return {this->x+q.x,this->y+q.y};
    }
    Vetor operator-(Vetor q) const{
        return {this->x-q.x,this->y-q.y};
    }
    Vetor operator*(double k) const{//Escalar * vetor
        return {this->x*k,this->y*k};
    }
    Vetor operator/(double k) const{//Vetor / escalar
        return {this->x/k,this->y/k};
    }
    double operator*(Vetor q) const{//Produto escalar
        return this->x*q.x+this->y*q.y;
    }
    double operator^(Vetor q) const{//Produto vetorial
        return this->x*q.y-q.x*this->y;
    }
    double operator~() const{//Modulo / distancia da origem
        return sqrt((*this)*(*this));
    }
    double dist(Vetor p,Vetor q) const{//Distancia de *this ao segmento pq
        Vetor r=*this;
        if((q-p)*(r-p)<=0||(p-q)*(r-q)<=0) return min(~(r-p),~(r-q));
        else return abs(((r-p)^(q-p))/(~(q-p)));
    }
    bool operator==(const Vetor p)const {
        return fabs(this->x - p.x) < EPS && fabs(this->y - p.y) < EPS;
    }
    bool operator!=(const Vetor p)const {
        return !(*this==p);
    }
    bool operator<(const Vetor p)const {
        if(this->x!=p.x) return this->x<p.x;
        return this->y<p.y;
    }
    static Vetor pivot;
    static bool cmpPolar(const Vetor& a, const Vetor& b){//Ordena por angulo polar a partir de pivot
        Vetor A = {a.x - pivot.x, a.y - pivot.y};
        Vetor B = {b.x - pivot.x, b.y - pivot.y};
        ll cross = A ^ B;
        if(cross == 0)
            return (A.x*A.x + A.y*A.y) < (B.x*B.x + B.y*B.y);
        return cross < 0; // <0 horario, >0 anti-horario
    }
};
Vetor Vetor::pivot = {0,0};

Vetor centroide(const vector<Vetor>& pts){
    ll sx=0, sy=0;
    for(auto &p: pts){ sx += p.x; sy += p.y; }
    return {sx/(ll)pts.size(), sy/(ll)pts.size()};
}

class Circle{
    public:
    Vetor o;
    double r;
    Circle(Vetor p,double raio): o(p),r(raio){ }
    Circle(double x,double y,double raio): o({x,y}),r(raio){ }
};

int ccw(Vetor a, Vetor b, Vetor c) {//1=esquerda, -1=direita, 0=colinear
    double cross = (b-a)^(c-a);
    if(fabs(cross) < EPS) return 0;
    return cross > 0 ? 1 : -1;
}
