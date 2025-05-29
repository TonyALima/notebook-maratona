
## Extra

### Utilidades

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
#define BitFlip(var,bit) var ^= (1<<bit)

// Maximo divisor comum (GCD) O(log10(min(a, b))
std::gcd(a, b);

//Ler linha com espaço
char str[500];
scanf("%[^\n]", str);
string s(str);
printf("%s",s.c_str());

//Transformar iterator em indice
distance(v.begin(),it);

```

<div style="page-break-after: always;"></div>

## Geometria

### Convex-Hull

```cpp
//Cria o menor poligono convexo que contem todos os pontos
//O(nlogn)
struct pt {
    double x, y;
    bool operator == (pt const& t) const {
        return x == t.x && y == t.y;
    }
};

int orientation(pt a, pt b, pt c) {
    double v = a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y);
    if (v < 0) return -1; // clockwise
    if (v > 0) return +1; // counter-clockwise
    return 0;
}

bool cw(pt a, pt b, pt c, bool include_collinear) {
    int o = orientation(a, b, c);
    return o < 0 || (include_collinear && o == 0);
}
bool collinear(pt a, pt b, pt c) { return orientation(a, b, c) == 0; }

void convex_hull(vector<pt>& a, bool include_collinear = false) {
    pt p0 = *min_element(a.begin(), a.end(), [](auto a, auto b) {
        return make_pair(a.y, a.x) < make_pair(b.y, b.x);
    });
    sort(a.begin(), a.end(), [&p0](const auto& a, const auto& b) {
        int o = orientation(p0, a, b);
        if (o == 0)
            return (p0.x-a.x)*(p0.x-a.x) + 
                    (p0.y-a.y)*(p0.y-a.y)
                    < (p0.x-b.x)*(p0.x-b.x) + 
                    (p0.y-b.y)*(p0.y-b.y);
        return o < 0;
    });
    if (include_collinear) {
        int i = (int)a.size()-1;
        while (i >= 0 && collinear(p0, a[i], a.back())) i--;
        reverse(a.begin()+i+1, a.end());
    }

    vector<pt> st;
    for (int i = 0; i < (int)a.size(); i++) {
        while (st.size() > 1 && 
                !cw(st[st.size()-2], st.back(), a[i], include_collinear))
            st.pop_back();
        st.push_back(a[i]);
    }

    if (include_collinear == false && st.size() == 2 && st[0] == st[1])
        st.pop_back();

    a = st;
}

```

<div style="page-break-after: always;"></div>

### Interseccao

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

## Grafos

### BFS

```cpp
//O(n^2) se mudar para lista de adjacencia vira O(n+a) sendo a o numero de arestas
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

### Bellman-Ford

```cpp
#define INF 0x3F3F3F3F
#define NMAX 100
// Caminho minimo com aresta negativa, caminho percorrido
//O(n^3)
int m[NMAX][NMAX], custo[NMAX], anterior[NMAX];

// Retorna true se existir ciclo negativo
bool bellmanFord(int s, int n){
    int i, j, k;
    for (i = 0; i < n; i++){
        custo[i] = INF;
        anterior[i] = -1;
    }
    custo[s] = 0;

    for (k = 0; k < n; k++)
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                if (m[i][j] != 0 && custo[i] != INF && custo[j] > custo[i] + m[i][j]) {
                    custo[j] = custo[i] + m[i][j];
                    anterior[j] = i;
                }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (m[i][j] != 0 && custo[i] != INF && custo[j] > custo[i] + m[i][j])
                return true; // Ciclo negativo
    
    return false;
}

```

<div style="page-break-after: always;"></div>

### DFS

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

### Dijkstra

```cpp
// Menor caminho
//O(e+nlogn) sendo e=arestas
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

### Edmonds-Karp

```cpp
//Usa BFS
//Necessario fonte e sumidouro
//O(ve^2) sendo v=vertices e=arestas
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

### Floyd-Warshall

```cpp
// Deteccao de ciclo negativo e caminho minimo para qualquer u, v
//O(n^3)
int m[][], custo[][];

bool floydWarshall(int n){
    int i, j, k;
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i == j) custo[i][j] = 0;
            else custo[i][j] = m[i][j];
        }
    }

    for (k = 0; k < n; k++)
        for (i = 0; i < n; i++){
            for (j = 0; j < n; j++){
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

### Kruskall

```cpp
// arvore geradora minima
// usa union find
//O(elogv)
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

### PrimDenso

```cpp
// árvore geradora mínima otimizada para grafos densos
// muito eficiente para a árvore geradora mínima dadas as posições cartesianas dos v vértices
// O(v²)
int n;
vector<vector<double>> adj; // MATRIZ de adjacência do grafo (nxn)
const double INF = 1000000000.0;

struct Edge
{
    double w = INF;
    int to = -1;
};

double prim()
{
    double total_weight = 0;
    vector<bool> selected(n, false);
    vector<Edge> min_e(n);
    min_e[0].w = 0;

    for (int i = 0; i < n; ++i)
    {
        int v = -1;
        for (int j = 0; j < n; ++j)
        {
            if (!selected[j] && (v == -1 || min_e[j].w < min_e[v].w))
                v = j;
        }

        /*if (min_e[v].w == INF) {
            cout << "No MST!" << endl;
            exit(0);
        }*/
        selected[v] = true;
        total_weight += min_e[v].w;
        /*if (min_e[v].to != -1)
            cout << v << " " << min_e[v].to << endl;*/

        for (int to = 0; to < n; ++to)
        {
            if (adj[v][to] < min_e[to].w)
                min_e[to] = {adj[v][to], v};
        }
    }

    return total_weight;
}
```

<div style="page-break-after: always;"></div>

### Tarjan

```cpp
#define NMAX 100
//O(v+e)
// componente fortemente conexo
vector<int> adj[NMAX];
stack<int> s;
int c, t, dis[NMAX], low[NMAX], ins[NMAX];
set<int> comp[NMAX]; set<int>::iterator sit;
// c eh o numero de componentes e comp[] contem cada componente

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
        c++;
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

### UnionFind

```cpp
#define MAXN 10000
//O(1) ná media, pior caso O(logn)
int parent[MAXN];
int w[MAXN];

int find_set(int v)
{
    if (v == parent[v])
        return v;
    return parent[v] = find_set(parent[v]);
}

void make_set(int v)
{
    parent[v] = v;
    w[v] = 1;
}

void union_sets(int a, int b)
{
    a = find_set(a);
    b = find_set(b);
    if (a != b)
    {
        if (w[a] < w[b])
            swap(a, b);
        parent[b] = a;
        w[a] += w[b];
    }
}

```

<div style="page-break-after: always;"></div>

## Matematica

### Busca-Ternaria

```cpp

const double EPS = 1e-6;
// Busca ternária para encontrar o ponto de máximo (ou mínimo, dependendo da função)
//O(logn) n = tamanho de busca
double ternary_search(std::function<double(double)> f, double left, double right) {
    while (right - left > EPS) {
        double m1 = left + (right - left) / 3.0;
        double m2 = right - (right - left) / 3.0;
        if (f(m1) < f(m2))  // Para máximo
            left = m1;
        else
            right = m2;
    }
    return (left + right) / 2.0;
}

/*Para numeros inteiros
int ternary_search(std::function<int(int)> f, int left, int right) {
    while (left<right) {
        int m1 = left + (right - left) / 3;
        int m2 = right - (right - left) / 3;
        if (f(m1) > f(m2))  // Para min
            left = m1+1;
        else
            right = m2-1;
    }
    return left;
}
*/

int main() {
    // Exemplo de função unimodal: f(x) = - (x - 2)^2 + 4, máximo em x = 2
    auto f = [](double x) {
        return -(x - 2) * (x - 2) + 4;
    };

    double a = 0.0, b = 4.0;//Limites superior e inferior

    // Encontra o ponto de máximo ou
    double xm = ternary_search(f, a, b);

    return 0;
}

```

<div style="page-break-after: always;"></div>

### Crammer

```cpp

//Usa determinante para resolver sistemas
//O(n^4)
double A[20][20],b[20],X[20];
//Calcula o determinante sem comprometer a matriz O(n^3)
double det(int n){
    double mat[n][n],res=1;
    int mxI;
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++) mat[i][j]=A[i][j];
    for(int i=0;i<n;i++){
        mxI=i;
        for(int j=i+1;j<n;j++){
            if(mat[j][i]>mat[mxI][i]) mxI=j;
        }
        if(mxI!=i){
            for(int j=0;j<n;j++){
                double aux=mat[mxI][j];
                mat[mxI][j]=mat[i][j];
                mat[i][j]=aux;
            }
            res=-res;
        }
        if(abs(mat[i][i])<1e-12) return NAN;
        for(int j=i+1;j<n;j++){
            double F=-mat[j][i]/mat[i][i];
            for(int k=0;k<n;k++) mat[j][k]+=F*mat[i][k];
        }
    }
    for(int i=0;i<n;i++) res*=mat[i][i];
    return res;
}

void solver(int n){
    double d=det(n);
    double aux[n];
    //if(abs(d)<1e-12) return;//determinante = 0
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            aux[j]=A[j][i];
            A[j][i]=b[j];
        }
        X[i]=det(n)/d;
        for(int j=0;j<n;j++) A[j][i]=aux[j];
    }
}
```

<div style="page-break-after: always;"></div>

### Crivo-de-Eratostenes

```cpp
// Cria um vetor de fatores primos / todos os primos ate N
//O(n)
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

### Crt

```cpp
struct Congruence {
    long long a, m;
};
//Calcular Chinese Remainder Theorem, usa inverso modular
long long mod_inv(long long a, long long m) {//O(m)
    if (a <= 1) return a;
    return m - (mod_inv(m % a, a) * (m / a) % m);
}

long long chinese_remainder_theorem(vector<Congruence> const& congruences) {//O(m*sizeof(congruences)
    long long M = 1;
    for (auto const& congruence : congruences) {
        M *= congruence.m;
    }

    long long solution = 0;
    for (auto const& congruence : congruences) {
        long long a_i = congruence.a;
        long long M_i = M / congruence.m;
        long long N_i = mod_inv(M_i, congruence.m);
        solution = (solution + a_i * M_i % M * N_i) % M;
    }
    return solution;
}

```

<div style="page-break-after: always;"></div>

### FFT

```cpp
using cd = complex<double>;
const double PI = acos(-1);
//Fast Fourier Transform O(nlogn)
void fft(vector<cd> & a, bool invert) {//Usado para multiplicar polinomios
    int n = a.size();

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (cd & x : a)
            x /= n;
    }
}
//Multiplicar polinomios e tambem somas possiveis de dois conjuntos de inteiros
vector<int> multiply(vector<int> const& a, vector<int> const& b) {
    vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size()) 
        n <<= 1;
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++)
        fa[i] *= fb[i];
    fft(fa, true);

    vector<int> result(n);
    for (int i = 0; i < n; i++)
        result[i] = round(fa[i].real());
    return result;
}

/*
int carry = 0;
    for (int i = 0; i < n; i++)
        result[i] += carry;
        carry = result[i] / 10;
        result[i] %= 10;
    }
*/

```

<div style="page-break-after: always;"></div>

### Fastpow

```cpp
//Potenciação rapida
//O(logn) n sendo power
long long fast_power(long long base, long long power) {
    long long result = 1;
    while(power > 0) {

        if(power % 2 == 1) { // Can also use (power & 1) to make code even faster
            result = (result*base) ;
        }
        base = (base * base);
        power = power / 2; // Can also use power >>= 1; to make code even faster
    }
    return result;
}

```

<div style="page-break-after: always;"></div>

### Gauss

```cpp
#define NMAX 20
//O(n^3)
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

## Programacao-Dinamica

### EditDistance

```cpp
//O(n*m)
int min(int a, int b, int c)
{
    if (a < b && a < c)
        return a;
    if (b < c)
        return b;
    return c;
}

int editDistance(string s1, string s2)
{
    int m = s1.size();
    int n = s2.size();
    int PD[m + 1][n + 1];
    // Inicialização da primeira coluna (E(i, 0))
    for (int i = 0; i <= m; i++)
        PD[i][0] = i;
    // Inicialização da primeira linha (E(0, j))
    for (int j = 0; j <= n; j++)
        PD[0][j] = j;
    // Preenchendo a matriz dp
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            if (s1[i - 1] == s2[j - 1])
            {
                PD[i][j] = PD[i - 1][j - 1]; // Sem custo extra
            }
            else
            {
                PD[i][j] = 1 + min(PD[i - 1][j], PD[i][j - 1], PD[i - 1][j - 1]);
            }
        }
    }
    return PD[m][n]; // Resposta final
}

```

<div style="page-break-after: always;"></div>

### Fenwick

```cpp
// range query soma
//O(logn)
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

### KnapSack

```cpp
// Encher a mochila com maior valor
//O(n*W)
int peso[], valor[];

// com repeticao
int knapSack(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int w = 0; w < W; w++){
        for (int i = 0; i < n; i++){
            if (peso[i] <= w)
                memo[w] = max(memo[w], memo[w-peso[i]] + valor[i]);
        }
    }
    return memo[W];
}

// sem repeticao
int knapSack(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int i = 0; i < n; i++){
        for (int w = W; w >=peso[i]; w--){
            memo[w] = max(memo[w], memo[w-peso[i]] + valor[i]);
        }
    }
    return memo[W];
}

// versao com matriz.

int mat[][];

int knapSack(int W, int n) {
    for (int i = 1; i <= n; ++i) {
        for (int w = 0; w <= W; ++w) {
            if (peso[i-1] <= w) {
                mat[i][w] = max(mat[i-1][w], mat[i-1][w - peso[i-1]] + valor[i-1]);
            } else {
                mat[i][w] = mat[i-1][w];
            }
        }
    }
    return mat[n][W];
}

vector<int> escolhidos(int W, int n){
    int w = W;
    vector<int> itens = vector<int>();
    while(n > 0 && w > 0){
        if (mat[n][w] != mat[n-1][w]) {
            itens.push_back(n-1);
            w -= peso[n-1];
        }
        n--;
    }
    return itens;
}

```

<div style="page-break-after: always;"></div>

### Lis-On-Tree

```cpp
//O(nlogn)
#define MAXN 100100
const int INF = 1e9;

int a[MAXN];
vector<int> tree[MAXN];
int r[MAXN];
int vis[MAXN];
vector<int> dp(MAXN, INF);
void dfs(int v)
{
    vis[v] = 1;
    int l = upper_bound(dp.begin(), dp.end(), a[v]) - dp.begin(); // índice do primeiro elemento em dp maior que a[i]
    int ant = dp[l];
    if (dp[l - 1] < a[v] && a[v] < dp[l])
    {
        dp[l] = a[v];
    }
    for (auto u : tree[v])
    {
        if (!vis[u])
        {
            dfs(u);
        }
    }
    r[v] = lower_bound(dp.begin(), dp.end(), INF) - dp.begin() - 1; // comprimento da LIS da raiz ao vértice v
    dp[l] = ant;
}

```

<div style="page-break-after: always;"></div>

### Lis

```cpp
//Calcula a maior subsequencia crescente
//O(nlogn)
int LIS(vector<int>& arr)
{
    // Binary search approach
    int n = arr.size();
    vector<int> ans;
    ans.push_back(arr[0]);
    for(int i=1;i<arr.size();i++){
        if (arr[i] > ans.back()) {
            // Se numero foi maior bota ele na sequencia
            ans.push_back(arr[i]);
        }
        else {
            //Senao troca o menor numero maior q ele por ele
            int low = lower_bound(ans.begin(), ans.end(),
                                  arr[i])
                      - ans.begin();
            ans[low] = arr[i];
        }
    }
    return ans.size();
}

int main(){
    int n;
    vector<int> v;
    cin>>n;
    while(n--){
        int aux;
        cin>>aux;
        v.push_back(aux);
    }
    cout<<LIS(v)<<endl;
    return 0;
}

```

<div style="page-break-after: always;"></div>

### MinimoDireita

```cpp
//O(n)
long long int vet[412345];
int mem[412345];
void dp(int n){
    stack<int> s;
    for(int i=n-1;i>=0;i--){
        while(!s.empty()&&vet[i]<=vet[s.top()]) s.pop();
        if(s.empty()) mem[i]=-1;
        else mem[i]=s.top();
        s.push(i);
    }
}

int main(){
    int n;
    long long int k;
    cin>>n>>k;
    for(int i=0;i<n;i++){
        cin>>vet[i];
        vet[i+n]=vet[i]-k*(n+i);
        vet[i]-=k*i;
    }
    dp(2*n);
    cout<<mem[0]%n+1;
    for(int i=1;i<n;i++){
        cout<<" "<<mem[i]%n+1;
    }
    cout<<endl;
    return 0;
}

```

<div style="page-break-after: always;"></div>

### Segment-Tree-Lazy-Propagation

```cpp
// range query/update, soma
//O(n)
const int maxn = 1000100;
int vet[maxn], tree[4 * maxn], lz[4 * maxn];

void build(int node, int l, int r)
{
    if (l == r)
    {
        tree[node] = vet[l];
        return;
    }
    int mid = (l + r) / 2;
    build(2 * node + 1, l, mid);
    build(2 * node + 2, mid + 1, r);
    tree[node] = tree[2 * node + 1] + tree[2 * node + 2];
}

void unlazy(int node, int tl, int tr)
{
    if (lz[node] == 0)
        return;
    tree[node] += (tr - tl + 1) * lz[node];
    if (tl != tr)
    {
        lz[2 * node + 1] += lz[node];
        lz[2 * node + 2] += lz[node];
    }
    lz[node] = 0;
}

void update(int node, int tl, int tr, int l, int r, int v)
{
    unlazy(node, tl, tr);
    if (tl > r || tr < l)
        return;
    if (tl >= l && tr <= r)
    {
        lz[node] += v;
        unlazy(node, tl, tr);
        return;
    }
    int mid = (tl + tr) / 2;
    update(2 * node + 1, tl, mid, l, r, v);
    update(2 * node + 2, mid + 1, tr, l, r, v);
    tree[node] = tree[2 * node + 1] + tree[2 * node + 2];
}

int query(int node, int tl, int tr, int l, int r)
{
    unlazy(node, tl, tr);
    if (r < tl || l > tr)
        return 0;
    if (l <= tl && r >= tr)
        return tree[node];
    int mid = (tl + tr) / 2;
    return query(2 * node + 1, tl, mid, l, r) +
           query(2 * node + 2, mid + 1, tr, l, r);
}

int main()
{
    int n = 5; // tamanho do vetor
    for (int i = 0; i < n; i++)
        vet[i] = 0;
    build(0, 0, n - 1);
    // Incrementa 3 no intervalo [1, 3]
    update(0, 0, n - 1, 1, 3, 3);
    // Soma total do intervalo [0, 4] -> deve ser 9 (0+3+3+3+0)
    cout << "Soma de [0, 4]: " << query(0, 0, n - 1, 0, 4) << '\n';
    return 0;
}

```

<div style="page-break-after: always;"></div>

### Segment-Tree

```cpp
// range query multiplicacao
//O(nlogn)
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

```

<div style="page-break-after: always;"></div>

## String

### Levenshtein

```cpp
//O(n*m)
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

### Trie

```cpp
//O(W*l) sendo W=numero de palavra l=tamanho da palavra
//Cria um dicionario de prefixos
const int K = 26;
struct Vertex {
    int next[K];
    int output = 0;
    bool eliminado = false;
    Vertex() {
        fill(begin(next), end(next), -1);
    }
};
vector<Vertex> trie(1);
//Adiciona o prefixo e da update no output
int add_string(string const& s) {
    int v = 0;
    for (char ch : s) {
        int c = ch - 'a';
        if(trie[v].eliminado) return 0;
        if (trie[v].next[c] == -1) {
            trie[v].next[c] = trie.size();
            trie.emplace_back();
        }
        v = trie[v].next[c];
    }
    if(trie[v].eliminado) return 0;
    trie[v].output++;
    return 1;
}
//Conta quantos prefixos terminam dps do ponto
int conta_out(int v){
    int res=0;
    if(trie[v].eliminado) return 0;
    for(int i=0;i<K;i++){
        if(trie[v].next[i]!=-1){
            res+=conta_out(trie[v].next[i]);
        }
    }
    res+=trie[v].output;
    return res;
}
//Apaga o prefixo e conta os eliminados
int remove_string(string const& s){
    int v = 0;
    for (char ch : s) {
        int c = ch - 'a';
        if(trie[v].eliminado) return 0;
        if (trie[v].next[c] == -1) {
            trie[v].next[c] = trie.size();
            trie.emplace_back();
        }
        v = trie[v].next[c];
    }
    int res=conta_out(v);
    trie[v].eliminado = true;
    return res;
}


int main(){
    int n,res=0;
    cin>>n;
    while(n--){
        int tipo;
        string s;
        cin>>tipo;
        cin>>s;
        if(tipo==2) res+=add_string(s);
        else res-=remove_string(s);
        cout<<res<<endl;
    }
    return 0;
}

```

<div style="page-break-after: always;"></div>

## Formulas Uteis

Soma de elementos em uma PA.

$S_n = E_1 + E_{n-1} \cdot \dfrac{n}{2}$

Soma de elementos em uma PG.

$Sn = a_1 \cdot \dfrac{(q^n - 1)}{q - 1}$
