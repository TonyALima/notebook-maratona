#include <bits/stdc++.h>
using namespace std;

// Requer: Template (Vetor, EPS, Circle, Reta)
// Obs: todas as funcoes aqui envolvem sqrt (interseccoes, tangentes, circuncentro),
// entao usam long double mesmo se Vetor virar long long int.
// Excecao: pointInCircle pode virar 100% inteira comparando distancia^2 com r^2
// (troque "~(p-c.o) <= c.r+EPS" por "(p-c.o)*(p-c.o) <= c.r*c.r" com Vetor/r em long long).
class Circle{
    public:
    Vetor o;
    long double r;//Pode ser trocado por long long int junto com Vetor::x,y
    Circle(Vetor p,long double raio): o(p),r(raio){ }
    Circle(long double x,long double y,long double raio): o({x,y}),r(raio){ }
    bool inCircle(Vetor p){//true se p esta dentro ou na borda
        return ~(p-o) <= r+EPS;
    }
    long double area(){
        return PI*r*r;
    }
    long double perimeter(){
        return 2*PI*r;
    }
    
    long double arco(long double angle){
        return fabsl(angle*r);
    }
    long double arco(Vetor a,Vetor b){
        long double cross = (b-o)^(a-o);
        long double dot = (a-o)*(b-o);
        return this->arco(atan2l(fabsl(cross),dot));
    }
};

bool pointInCircle(Circle c, Vetor p){//true se p esta dentro ou na borda de c
    return ~(p-c.o) <= c.r+EPS;
}

vector<Vetor> lineCircleIntersection(Reta r, Circle c){//0, 1 ou 2 pontos de interseccao
    Vetor n=~r, dir=!r;
    long double d=r.posicao(c.o)/(~n); // distancia (com sinal) do centro a reta
    if(fabsl(d)>c.r+EPS) return {};
    Vetor proj=c.o-(n/(~n))*d; // pe da perpendicular do centro na reta
    long double h=sqrtl(max((long double)0.0,c.r*c.r-d*d));
    Vetor u=dir/(~dir);
    if(h<EPS) return {proj};
    return {proj+u*h, proj-u*h};
}

vector<Vetor> circleCircleIntersection(Circle a, Circle b){//0, 1, 2 pontos (ou infinitos se coincidentes, retorna vazio)
    Vetor d=b.o-a.o;
    long double dist=~d;
    if(dist<EPS) return {}; // concentricos: infinitas ou nenhuma interseccao
    if(dist>a.r+b.r+EPS || dist<fabsl(a.r-b.r)-EPS) return {};
    long double x=(dist*dist-b.r*b.r+a.r*a.r)/(2*dist); // distancia do centro a ate o ponto medio da corda
    long double h2=a.r*a.r-x*x;
    long double h=sqrtl(max((long double)0.0,h2));
    Vetor u=d/dist, n={-u.y,u.x};
    Vetor p=a.o+u*x;
    if(h<EPS) return {p};
    return {p+n*h, p-n*h};
}

Circle circumcircle(Vetor a, Vetor b, Vetor c){//circulo que passa por a, b e c; UB se a,b,c colineares (d==0)
    Vetor ab=b-a, ac=c-a;
    long double d=2*(ab^ac); // d==0 se colineares: verificar antes de chamar
    long double ab2=ab*ab, ac2=ac*ac;
    Vetor o={a.x+(ac.y*ab2-ab.y*ac2)/d, a.y+(ab.x*ac2-ac.x*ab2)/d};
    return Circle(o, ~(o-a));
}

Circle incircle(Vetor a, Vetor b, Vetor c){
    long double r=fabsl((a-b)^(b-c))/(~(a-b)+~(b-c)+~(c-a));
    Vetor o=((c*~(a-b))+(a*~(b-c))+(b*~(c-a)))/(~(a-b)+~(b-c)+~(c-a));
    return Circle(o,r);
}
