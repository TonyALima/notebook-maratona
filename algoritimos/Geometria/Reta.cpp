#include <bits/stdc++.h>
using namespace std;

// Requer: Template (Vetor, EPS)
// Reta: ~=normal unitario !=direcao unitario posicao(p) dist(p) ||=paralelo ==igual ^=interseccao
// a,b,c sao combinacoes lineares de coordenadas: se Vetor::x,y forem ll, troque a,b,c para ll
// (posicao() fica exata). dist() e operator^ envolvem divisao/sqrt e devem seguir em long double.
class Reta{//ax+by+c=0 | construida a partir de dois pontos
    public:
    long double a,b,c;
    Reta(Vetor s,Vetor e){
        this->a=s.y-e.y;
        this->b=e.x-s.x;
        this->c=s^e;
    }
    Vetor operator~() const{//Vetor normal
        Vetor n={this->a,this->b};
        return n;
    }
    Vetor operator!() const{//Vetor direcao
        Vetor d={-this->b,this->a};
        return d;
    }
    long double posicao(Vetor p){//=0 na reta, >0 ou <0 nos lados
        return this->a*p.x+this->b*p.y+c;
    }
    long double dist(Vetor p){
        return fabsl(this->posicao(p))/(~(~(*this)));
    }
    bool operator||(const Reta r)const {//true se paralelas ou iguais
        return (fabsl(this->a*r.b - this->b*r.a) < EPS);
    }
    bool operator==(const Reta& l) const {
        return fabsl(this->a*l.b - this->b*l.a) < EPS &&
            fabsl(this->a*l.c - this->c*l.a) < EPS &&
            fabsl(this->b*l.c - this->c*l.b) < EPS;
    }
    bool intersect(Vetor a,Vetor b){//Reta intersecta segmento
        long double pa=this->posicao(a),pb=this->posicao(b);
        if((pa>0&&pb<0)||(pa<0&&pb>0)||(fabs(pa)<EPS||fabs(pb)<EPS)) return true;
        return false;
    }
    long double dist(Vetor p,Vetor q){//Distancia da reta ao segmento pq
        if(this->intersect(p,q)) return 0;
        return min(this->dist(p),this->dist(q));
    }
    long double dist(Reta r){//Distancia entre retas
        if(!(*this==r) && (*this||r)){
            long double d2=this->a*this->a+this->b*this->b;
            Vetor p={-this->a*this->c/d2,-this->b*this->c/d2};//ponto qualquer desta reta
            return r.dist(p);
        }
        return 0;
    }
    Vetor operator^(Reta r){//Ponto de interseccao (verificar || antes)
        long double det=this->a*r.b-this->b*r.a;
        Vetor ponto={-(this->c*r.b-this->b*r.c)/det,-(this->a*r.c-this->c*r.a)/det};
        return ponto;
    }
};
