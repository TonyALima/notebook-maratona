#include <bits/stdc++.h>
using namespace std;

// Requer: Template (Vetor, EPS, ccw), Reta e Interseccao (onSegment)
// Raio: semirreta a partir de o na direcao d (nao precisa ser unitaria)
// Se Vetor::x,y forem ll, troque "fabsl(v^d)>EPS" por "(v^d)!=0" (comparacao exata).
class Raio{
    public:
    Vetor o,d;
    Raio(Vetor origem, Vetor ponto): o(origem), d(ponto-origem) { }
    long double dist(Vetor p) const{//distancia de p ate a semirreta
        Vetor v=p-o;
        if((v*d)<=0) return ~v; // projecao cai antes de o: ponto mais proximo eh o
        Reta r(o,o+d);
        return fabsl(r.dist(p));
    }
    long double dist(Raio smr) const{
        Reta r1(o,o+d),r2(smr.o,smr.o+smr.d);
        if(r1==r2){
            if((d*smr.d)>0) return 0; // mesma direcao: sempre se sobrepoem
            if((smr.o-o)*d>=-EPS) return 0; // direcoes opostas mas se alcancam
            return ~(o-smr.o); // direcoes opostas, gap entre as origens
        }
        if(r1||r2) return min(this->dist(smr.o),smr.dist(o));
        Vetor p=r1^r2;
        if((p-o)*d>=-EPS && (p-smr.o)*smr.d>=-EPS) return 0;
        return min(this->dist(smr.o),smr.dist(o));
    }
    long double dist(Reta r) const{
        Reta r1(o,o+d);
        if(r1||r) return r.dist(o);
        Vetor p=r^r1;
        if((p-o)*d<=0) return  r.dist(o);
        else return 0;
    }
    bool intersect(Vetor p,Vetor q) const{//semirreta intersecta segmento pq
        Reta r1(o,o+d),r2(p,q);
        if(r1||r2){
            if(!(r1==r2)) return false;
            return (p-o)*d>=-EPS || (q-o)*d>=-EPS || onSegment(p,q,o);
        }
        Vetor x=r1^r2;
        if((x-o)*d<-EPS) return false;
        return onSegment(p,q,x);
    }
    long double dist(Vetor p,Vetor q) const{//Distancia da semirreta ao segmento pq
        if(this->intersect(p,q)) return 0;
        return min({this->dist(p),this->dist(q),o.dist(p,q)});
    }
};
