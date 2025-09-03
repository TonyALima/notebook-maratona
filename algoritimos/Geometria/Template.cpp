#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;
const double EPS = 1e-9;
const double PI = acos(-1);

class Vetor{
    public: 
    double x,y;//Pode ser trocado por long long int dependendo da questao
    Vetor operator+(Vetor q) const{
        return {this->x+q.x,this->y+q.y};
    }
    Vetor operator-(Vetor q) const{
        return {this->x-q.x,this->y-q.y};
    }
    Vetor operator*(double k) const{//Escalar vezes vetor
        return {this->x*k,this->y*k};
    }
    Vetor operator/(double k) const{//Escalar vezes vetor
        return {this->x/k,this->y/k};
    }
    double operator*(Vetor q) const{//Vetor escalar vetor
        return this->x*q.x+this->y*q.y;
    }
    double operator^(Vetor q) const{//Produto vetorial
        return this->x*q.y-q.x*this->y;
    }
    double operator~() const{//Modulo do vetor/Distancia
        return sqrt((*this)*(*this));
    }
    double dist(Vetor p,Vetor q) const{
        Vetor r=*this;
        if((q-p)*(r-p)<=0||(p-q)*(r-q)<=0) return min(~(r-p),~(r-q));
        else return abs(((r-p)^(q-p))/(~(q-p)));
    }
    bool operator==(const Vetor p)const {
        return fabs(this->x - p.x) < EPS && fabs(this->y - p.y) < EPS;
    }
    bool operator!=(const Vetor p)const {
        if(*this==p) return false;
        return true;
    }
    bool operator<(const Vetor p)const {
        if(this->x!=p.x) return this->x<p.x;
        return this->y<p.y;
    }
    static Vetor pivot; 

    static bool cmpPolar(const Vetor& a, const Vetor& b){//Usado no sort pra ordenar no sentido horario
        Vetor A = {a.x - pivot.x, a.y - pivot.y};
        Vetor B = {b.x - pivot.x, b.y - pivot.y};
        ll cross = A ^ B;
        if(cross == 0) {
            // se colineares, ordenar pelo mais próximo ao pivô
            return (A.x*A.x + A.y*A.y) < (B.x*B.x + B.y*B.y);
        }
        return cross < 0; // "<0" = sentido horário, ">0" = anti-horário
    }
};
Vetor Vetor::pivot = {0,0};//Vai receber o valor da funcao centroide
Vetor centroide(const vector<Vetor>& pts){//Calcula o pivot pra ordenar
    ll sx=0, sy=0;
    for(auto &p: pts){
        sx += p.x;
        sy += p.y;
    }
    return {sx/(ll)pts.size(), sy/(ll)pts.size()};
}

class Reta{//Caso seja passado pontos para descobrir a reta, resolver o sistema
    public:
    double a,b,c;//ax+by+c
    Reta(Vetor s,Vetor e){
        this->a=s.y-e.y;
        this->b=e.x-s.x;
        this->c=s^e;
    }
    Vetor operator~() const{//Vetor normal a reta
        Vetor n={this->a,this->b};
        return n/(~n);
    }
    Vetor operator!() const{//Vetor direcao da reta
        Vetor d={-this->b,this->a};
        return d/(~d);
    }
    double posicao(Vetor p){//Na reta =0, de um lado >0 no outro <0
        return this->a*p.x+this->b*p.y+c;
    }
    double dist(Vetor p){
        return (this->posicao(p))/(~(~(*this)));
    }
    bool operator&&(const Reta r)const {//Returna  true se for paralelo
        return (fabs(this->a*r.b - this->b*r.a) < EPS);
    }
    bool operator==(const Reta& l) const {
        return fabs(this->a*l.b - this->b*l.a) < EPS &&
            fabs(this->a*l.c - this->c*l.a) < EPS &&
            fabs(this->b*l.c - this->c*l.b) < EPS;
    }
    Vetor operator^(Reta r){//Retorna ponto de interseccao
        double det=this->a*r.b-this->b*r.a;//Teoricamente !=0 pois foi usado & antes
        Vetor ponto={-(this->c*r.b-this->b*r.c)/det,-(this->a*r.c-this->c*r.a)/det};
        return ponto;
    }
};

class Circle{
    public:
    Vetor o;
    double r;
    Circle(Vetor p,double raio): o(p),r(raio){ }
    Circle(double x,double y,double raio): o({x,y}),r(raio){ }
};

int ccw(Vetor a, Vetor b, Vetor c) {
    double cross = (b-a)^(c-a);
    if(fabs(cross) < EPS) return 0; // colineares
    return cross > 0 ? 1 : -1; // 1 = esquerda, -1 = direita | c da reta ab
}

bool onSegment(Vetor a, Vetor b, Vetor p) {//p no segmento ab
    return ccw(a,b,p)==0 && 
       min(a.x,b.x)-EPS <= p.x && p.x <= max(a.x,b.x)+EPS &&
       min(a.y,b.y)-EPS <= p.y && p.y <= max(a.y,b.y)+EPS;
}

bool intersect(Vetor a, Vetor b, Vetor c, Vetor d) {//segmento ab e cd intersectam
    int d1 = ccw(a,b,c), d2 = ccw(a,b,d);
    int d3 = ccw(c,d,a), d4 = ccw(c,d,b);
    if(d1*d2<0 && d3*d4<0) return true;
    if(d1==0 && onSegment(a,b,c)) return true;
    if(d2==0 && onSegment(a,b,d)) return true;
    if(d3==0 && onSegment(c,d,a)) return true;
    if(d4==0 && onSegment(c,d,b)) return true;
    return false;
}

double area(vector<Vetor>& poly){//Precisa estar ordenado no sentido horario ou contrario e ser convexo
    double s=0;
    int n=poly.size();
    for(int i=0;i<n;i++)
        s += poly[i] ^ poly[(i+1)%n];
    return fabs(s)/2;
}

int main(){

    return 0;
}
