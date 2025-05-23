#include <bits/stdc++.h>
using namespace std;
// range query multiplicacao

int vet[];
struct tree{
    tree *esq, *dir;
    int from, to, valor;
    tree(int _from, int _to):from(_from), to(_to), dir(NULL), esq(NULL), valor(1){}
};

tree * build(int e, int d){
    if (e > d) return NULL;
    tree *res = new tree(e,d);
    if (e == d) res->valor = vet[e];
    else{
        int m = (e+d)/2;
        res->esq = build(e, m); res->dir = build(m+1, d);
        if (res->esq != NULL) res->valor *= res->esq->valor;
        if (res->dir != NULL) res->valor *= res->dir->valor;
    }
    return res;
}

int query(tree *arv, int e, int d){
    if (arv == NULL) return 1;
    if (e <= arv->from && arv->to <= d) return arv->valor;
    if (e > arv->to) return 1; if (d < arv->from) return 1;
    return query(arv->esq, e, d)*query(arv->dir, e, d);
}

int update(tree *arv, int i, int valor){
    if (arv == NULL) return 1;
    if (arv->to < i || arv->from > i) return arv->valor;
    if (arv->from == arv->to && arv->from == i) arv->valor = valor;
    else arv->valor = (update(arv->esq, i, valor)*update(arv->dir, i, valor));
    return arv->valor;
}
