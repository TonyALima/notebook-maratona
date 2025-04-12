# Notebook de Maratona

## Tarjan

```cpp
// componente fortemente conexo

vector<int> adj[];
stack<int> s;
int c, t, dis[], low[], ins[];
set<int> comp[]; set<int>::iterator sit;

void DFS(int u){
    dis[u] = low[u] = t++;
    s.push(u); ins[u] = 1;
    vector<int>::iterator it;
    for(it = adj[u].begin(); it != adj[u].end(); it++){
        int v = *it;
        if(dis[v] == -1){ DFS(v); low[u] = min(low[u], low[v]);}
        else if(ins[v] == 1) low[u] = min(low[u], dis[v]);
    }
    if (low[u] == dis[u]){
        int w = 0;
        while (s.top() != u){
            w = s.top();
            comp[c].insert(w); ins[w] = 0; s.pop();
        }
        w = s.top();
        comp[c].insert(w); ins[w] = 0; s.pop();
    } 
}

void tarjan(int n){
    t = c = 0;
    for (int i = 0; i < n; i++){
        dis[i] = low[i] = -1; ins[i] = 0; comp[i].clear();
    }
    for(int i = 0; i < n; i++) if (dis[i] == -1) DFS(i);
}
```

## Edmonds-Karp

```cpp
//Usa BFS
//Necessario fonte e sumidouro

int fluxoMaximo(int ini, int fim, int tam){
    int u, v;
    int fluxo = 0;
    int bot;
    int rgrafo[tam][tam];
    for(u = 0; u < tam; u++)
        for(v = 0; v < tam; v++) rgrafo[u][v]=mat[u][v];
    
    while (BFS(ini, fim, tam)){
        bot = INF;
        for (v = fim; v != ini; v = anterior[v]){
            u = anterior[v];
            bot = min(bot, rgrafo[u][v]);
        }
        for (v = fim; v != ini; v = anterior[v]){
            u = anterior[v];
            rgrafo[u][v] -= bot; rgrafo[v][u] += bot; 
        }
        fluxo += bot;
    }
    return fluxo;
}
```

## BFS

```cpp
#define INF 0x3F3F3F3F

int mat[][];
int vis[], anterior[];

int BFS(int ini, int fim, int tam){ // 1 se tiver caminho, 0 caso nao
    int i;
    queue<int> fila;
    for (i = 0; i < tam; i++){ vis[i] = 0; }
    fila.push(ini); vis[ini] = 1; anterior[ini] = -1;
    while (!fila.empty()){
        for (i = 0; i < tam; i++){
            if (!vis[i] && mat[fila.front()][i]){
                if (i == fim) anterior[i] = fila.front(); return 1;
                fila.push(i); anterior[i] = fila.front(); vis[i] = 1;
            }
        }
        fila.pop();
    }
    return 0;
}
```

## Kruskall

```cpp
// arvore geradora minima
// usa union find

int mat[][];
int pai[], rnk[];

int find(int u){
    while (u != pai[u]) u = pai[u];
    return pai[u];
}

int merge(int u, int v){
    u = find(u); v = find(v);
    if(rnk[u] > rnk[v]) pai[v] = u; else pai[u] = v;
    if(rnk[u] == rnk[v]) rnk[v]++;
}

int kruskall(int tam){
    vector<tuple<int, int, int>> vec;
    vector<tuple<int, int, int>>::iterator it;
    int i, j, res = 0;
    for (i = 0; i < tam; i++){
        pai[i] = i; rnk[i] = 0;
        for (j = 0; j < tam; j++) vec.push_back(make_tuple(mat[i][j], i, j));
    }
    sort(vec.begin(), vec.end());
    for(it = vec.begin(); it != vec.end(); it++)
        if (find(get<1>(*it)) != find(get<2>(*it))){
            res += get<0>(*it); merge(get<1>(*it), get<1>(*it));
        }
    return res;
}
```

## Interseccao

```cpp

typedef struct{double x, y;}Ponto;

// Função para calcular o produto vetorial
double crossProduct(Ponto a, Ponto b) {
    return a.x * b.y - a.y * b.x;
}

bool doIntersect(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    double d1 = crossProduct({p3.x - p1.x, p3.y - p1.y}, {p2.x - p1.x, p2.y - p1.y});
    double d2 = crossProduct({p4.x - p1.x, p4.y - p1.y}, {p2.x - p1.x, p2.y - p1.y});
    double d3 = crossProduct({p1.x - p3.x, p1.y - p3.y}, {p4.x - p3.x, p4.y - p3.y});
    double d4 = crossProduct({p2.x - p3.x, p2.y - p3.y}, {p4.x - p3.x, p4.y - p3.y});

    return (d1 * d2 < 0) && (d3 * d4 < 0);
}
```

## Gauss

```cpp
#define NMAX 20

// matriz ampliada de coeficientes
double mat[NMAX][NMAX+1];
double sol[NMAX];

void elimination(int n){
    double factor;
    for (int i = 0; i < n - 1; i++){
        for (int j = i+1; j < n; j++){
            factor = mat[j][i]/mat[i][i];
            for (int k = 0; k < n + 1; k++){
                mat[j][k] = mat[j][k] - factor*mat[i][k];
            }
        }
    }
}

void solve(int n){
    memset(sol, 0, sizeof(sol));
    double s = 0;
    for (int l = n - 1; l >= 0; l--) {
        s = 0; 
        for (int j = l; j < n; j++) {
            s += sol[j]*mat[l][j];
        }
        sol[l] = (mat[l][n] - s) / mat[l][l];
    }
}
```

## DFS

```cpp
// deteccao de ciclos

int mat[][];
int vis[], rec[]; // inicializa com 0

int DFS(int v, int tam){
    if(!vis[v]){
        vis[v] = 1;
        rec[v] = 1;
        for (int i = 0; i < tam; i++){
            if (mat[v][i]){
                if(!vis[i] && DFS(i, tam)) return 1;
                else if(rec[i]) return 1;
            }
        }
    }
    rec[v] = 0;
    return 0;
}
```

## Fenwick

```cpp
// range query soma

int BIT[];

void updateBIT(int tam, int index, int valor){
    index++;
    while (index <= tam){
        BIT[index] += valor;
        index += index & (-index);
    }
}

void buildBIT(int  *vet, int tam){
    int i;
    memset(BIT, 0, sizeof(BIT));
    for (i == 0; i < tam; i++) updateBIT(tam, i, vet[i]);
}

int queryBIT(int index){
    int soma = 0;
    while(index > 0){
        soma += BIT[index];
        index -= index & (-index);
    }
    return soma;
}
```

## Segment-Tree

```cpp
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
    if (e <= arv->from && arv->to >= d) return arv->valor;
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
```

## Utilidades

```cpp

// Precisao 2 casas decimais de float para impressao
cout << fixed << setprecision(2);

// Criar pair
make_pair(1, 2);

// Criar tupla
make_tuple(1, 2, 3);

// Pegar elemento i da tupla
get<i>(t);
```
