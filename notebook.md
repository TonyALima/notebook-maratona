# Notebook de Maratona

## Levenshtein

```cpp

// numero minimo de operacoes para transformar uma string em outra.
int levenshtein(const string& a, const string& b) {
    int m = a.size(), n = b.size();
    vector<vector<int>> dp(m+1, vector<int>(n+1));

    // Inicializa bordas
    for (int i = 0; i <= m; ++i) dp[i][0] = i;
    for (int j = 0; j <= n; ++j) dp[0][j] = j;

    // Preenche a matriz
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int cost = (a[i-1] == b[j-1]) ? 0 : 1;
            dp[i][j] = min({
                dp[i-1][j] + 1,     // deletar
                dp[i][j-1] + 1,     // inserir
                dp[i-1][j-1] + cost // substituir
            });
        }
    }

    return dp[m][n];
}

```

<div style="page-break-after: always;"></div>

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

<div style="page-break-after: always;"></div>

## Edmonds-Karp

```cpp
//Usa BFS
//Necessario fonte e sumidouro
int rgrafo[][]; // alterar BFS para usar grafo residual.
int fluxoMaximo(int ini, int fim, int tam){
    int u, v;
    int fluxo = 0;
    int bot;
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

<div style="page-break-after: always;"></div>

## BFS

```cpp

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
                if (i == fim) {anterior[i] = fila.front(); return 1;}
                fila.push(i); anterior[i] = fila.front(); vis[i] = 1;
            }
        }
        fila.pop();
    }
    return 0;
}

```

<div style="page-break-after: always;"></div>

## Crivo-de-Eratostenes

```cpp
// Cria um vetor de fatores primos / todos os primos ate N

int fprimos[]; // inicializa com 0

void crivo(int n){
    fprimos[0] = fprimos[1] = 1;
    for(int i = 2; i < n; i++){
        if(fprimos[i] == 0){
            fprimos[i] = i;
        }
        if (i * i < n){
            for (int j = i*i; j < n; j+=i){
                fprimos[j] = i;
            }
        }
    }
}

```

<div style="page-break-after: always;"></div>

## Dijkstra

```cpp
// Menor caminho

typedef pair<int, int> ii;
typedef vector<ii> vii;
typedef vector<int> vi;

vii adj[];
vi custo(SIZE, INF);

void dijkstra(int s){
    custo[s] = 0;
    priority_queue<ii, vii, greater<ii>> heap;
    heap.push(ii(0, s));
    while(!heap.empty()){
        ii menor = heap.top(); heap.pop();
        int d = menor.first, u = menor.second;
        if (d > custo[u]) continue;
        for (int i = 0; i < adj[u].size(); i++){
            ii v = adj[u][i];
            if(custo[u]+v.second < custo[v.first]){
                custo[v.first] = custo[u] + v.second;
                heap.push(ii(custo[v.first], v.first));
            }
        }
    }
}

```

<div style="page-break-after: always;"></div>

## Bellman-Ford

```cpp
// Caminho minimo com aresta negativa, caminho percorrido

int m[][], custo[], anterior[];

void bellmanFord(int s, int n){
    int i, j, k;
    for (i = 0; i < n; i++){
        custo[i] = INF;
        anterior[i] = -1;
    }
    custo[s] = 0;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < n; k++)
                if (k != j && m[j][k] != 0 && custo[k] > custo[j]+m[j][k]){
                    custo[k] = custo[j]+m[j][k];
                    anterior[k] = j;
                }
}

```

<div style="page-break-after: always;"></div>

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

<div style="page-break-after: always;"></div>

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

<div style="page-break-after: always;"></div>

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

<div style="page-break-after: always;"></div>

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

<div style="page-break-after: always;"></div>

## Floiyd-Warshall

```cpp
// Deteccao de ciclo negativo e caminho minimo para qualquer u, v

int m[][], custo[][];

bool floydWarshall(int n){
    int i, j, k;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            custo[i][j] = m[i][j];

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++){
            for (k = 0; k < n; k++){
                if (custo[i][k] != INF && custo[k][j] != INF &&
                    custo[i][j] > custo[i][k] + custo[k][j])
                    custo[i][j] = custo[i][k] + custo[k][j];
            }
            if (custo[j][j]< 0)
                return true;
        }
    return false;
}

```

<div style="page-break-after: always;"></div>

## KnapSack

```cpp
// Encher a mochila com maior valor

pair<int, int> vet[]; // peso, valor

// sem repeticao
int knapSack(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int i = 0; i < n; i++){
        for (int w = W; w >=vet[i].first; w--){
            memo[w] = max(memo[w], memo[w-vet[i].first] + vet[i].second);
        }
    }
    return memo[W];
}

// com repeticao
int knapSack(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int w = 0; w < W; w++){
        for (int i = 0; i < n; i++){
            if (vet[i].first <= w)
                memo[w] = max(memo[w], memo[w-vet[i].first] + vet[i].second);
        }
    }
    return memo[W];
}

```

<div style="page-break-after: always;"></div>

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

<div style="page-break-after: always;"></div>

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

<div style="page-break-after: always;"></div>

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

// infinito
#define INF 0x3F3F3F3F

// Operacoes BitWise 
#define BitTest(var, bit) var & (1 << bit)
#define BitSet(var, bit) var |= (1 << bit)
#define BitClear(var, bit) var &= ~(1 << bit)

```

<div style="page-break-after: always;"></div>

## Formulas Uteis

Soma de elementos em uma PA.

$S_n = E_1 + E_{n-1} \cdot \dfrac{n}{2}$

Soma de elementos em uma PG.

$Sn = a_1 \cdot \dfrac{(q^n - 1)}{q - 1}$
