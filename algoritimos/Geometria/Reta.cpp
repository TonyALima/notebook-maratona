#include <bits/stdc++.h>
using namespace std;

// Requer: Template (Vetor, EPS)
// Reta: ~=normal unitario !=direcao unitario posicao(p) dist(p) ||=paralelo ==igual ^=interseccao
class Reta{//ax+by+c=0 | construida a partir de dois pontos
    public:
    double a,b,c;
    Reta(Vetor s,Vetor e){
        this->a=s.y-e.y;
        this->b=e.x-s.x;
        this->c=s^e;
    }
    Vetor operator~() const{//Vetor normal unitario
        Vetor n={this->a,this->b};
        return n/(~n);
    }
    Vetor operator!() const{//Vetor direcao unitario
        Vetor d={-this->b,this->a};
        return d/(~d);
    }
    double posicao(Vetor p){//=0 na reta, >0 ou <0 nos lados
        return this->a*p.x+this->b*p.y+c;
    }
    double dist(Vetor p){
        return (this->posicao(p))/(~(~(*this)));
    }
    bool operator||(const Reta r)const {//true se paralelas
        return (fabs(this->a*r.b - this->b*r.a) < EPS);
    }
    bool operator==(const Reta& l) const {
        return fabs(this->a*l.b - this->b*l.a) < EPS &&
            fabs(this->a*l.c - this->c*l.a) < EPS &&
            fabs(this->b*l.c - this->c*l.b) < EPS;
    }
    Vetor operator^(Reta r){//Ponto de interseccao (verificar && antes)
        double det=this->a*r.b-this->b*r.a;
        Vetor ponto={-(this->c*r.b-this->b*r.c)/det,-(this->a*r.c-this->c*r.a)/det};
        return ponto;
    }
};
