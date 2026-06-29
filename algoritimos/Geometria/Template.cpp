#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const long double EPS = 1e-9;
const long double PI = acosl(-1);

// Vetor: +(add) -(sub) *(long double)=escala *(Vetor)=dot ^=cross ~=modulo dist(p,q)=dist ao segmento pq
// cmpPolar: usar com sort + Vetor::pivot para ordenar por angulo polar
// Para coordenadas inteiras: troque x,y (e Circle::r) para ll. operator*, operator^ e cmpPolar
// continuam exatos em ll; operator~, dist e o construtor de Reta/funcoes de Circulo que usam
// sqrt seguem precisando de long double (faca o calculo em ll e converta na hora do sqrt).
class Vetor{
    public:
    long double x,y;//Pode ser trocado por long long int dependendo da questao
    Vetor operator+(Vetor q) const{
        return {this->x+q.x,this->y+q.y};
    }
    Vetor operator-(Vetor q) const{
        return {this->x-q.x,this->y-q.y};
    }
    Vetor operator*(long double k) const{//Escalar * vetor
        return {this->x*k,this->y*k};
    }
    Vetor operator/(long double k) const{//Vetor / escalar
        return {this->x/k,this->y/k};
    }
    long double operator*(Vetor q) const{//Produto escalar
        return this->x*q.x+this->y*q.y;
    }
    long double operator^(Vetor q) const{//Produto vetorial
        return this->x*q.y-q.x*this->y;
    }
    long double operator~() const{//Modulo / distancia da origem
        return sqrtl((*this)*(*this));
    }
    long double dist(Vetor p,Vetor q) const{//Distancia de *this ao segmento pq
        Vetor r=*this;
        // dot<=0: projecao cai antes de p (ou apos q), entao o extremo mais proximo eh p (ou q)
        if((q-p)*(r-p)<=0||(p-q)*(r-q)<=0) return min(~(r-p),~(r-q));
        else return fabsl(((r-p)^(q-p))/(~(q-p)));
    }
    bool operator==(const Vetor p)const {
        return fabs(this->x - p.x) < EPS && fabs(this->y - p.y) < EPS;
    }
    bool operator!=(const Vetor p)const {
        return !(*this==p);
    }
    bool operator<(const Vetor p)const {
        if(fabsl(this->x - p.x) >= EPS) return this->x < p.x;
        return this->y < p.y;
    }
    static Vetor pivot;
    static bool cmpPolar(const Vetor& a, const Vetor& b){//Ordena por angulo polar a partir de pivot
        Vetor A = {a.x - pivot.x, a.y - pivot.y};
        Vetor B = {b.x - pivot.x, b.y - pivot.y};
        long double cross = A ^ B; // se x,y forem ll, troque para "ll cross = A ^ B" (operacao exata)
        if(cross == 0)
            return (A.x*A.x + A.y*A.y) < (B.x*B.x + B.y*B.y);
        return cross < 0; // <0 horario, >0 anti-horario
    }
};
Vetor Vetor::pivot = {0,0};

Vetor centroide(const vector<Vetor>& pts){
    long double sx=0, sy=0;
    for(auto &p: pts){ sx += p.x; sy += p.y; }
    return {sx/pts.size(), sy/pts.size()};
}

int ccw(Vetor a, Vetor b, Vetor c) {//1=esquerda, -1=direita, 0=colinear
    long double cross = (b-a)^(c-a); // se x,y forem ll, troque para "ll cross" e "if(cross==0) return 0;" (sem EPS)
    if(fabsl(cross) < EPS) return 0;
    return cross > 0 ? 1 : -1;
}
