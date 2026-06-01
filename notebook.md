
## Extra

### Utilidades

```cpp

//Entrada rapida no cin
ios::sync_with_stdio(false);
cin.tie(nullptr);

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
//Vai ser utilizado a classe Vetor
//da classe Vetor - ^ == < pivot cmpPolar
//ccw
bool cw(Vetor a, Vetor b, Vetor c,bool include_linear){
    int o=ccw(a,b,c);
    return o<0||(include_linear&&o==0);
}

void convex_hull(vector<Vetor>& pts, bool include_collinear=false) {//O(nlogn)
    if(pts.size() <= 1) return;

    // escolher pivô (menor y, e se empatar menor x)
    Vetor::pivot = *min_element(pts.begin(), pts.end(), [](const Vetor& a, const Vetor& b){
        return make_pair(a.y, a.x) < make_pair(b.y, b.x);
    });

    // ordenar pelo ângulo polar usando a função do template
    sort(pts.begin(), pts.end(), Vetor::cmpPolar);

    if (include_collinear) {
        int i = (int)pts.size()-1;
        while (i >= 0 && ccw(Vetor::pivot, pts[i], pts.back()) == 0) i--;
        reverse(pts.begin()+i+1, pts.end());
    }

    // construção da casca convexa
    vector<Vetor> st;
    for (int i=0; i<(int)pts.size(); i++) {
        while (st.size() > 1 && !cw(st[st.size()-2],st.back(),pts[i],include_collinear)) // usa cw para ficar mais legivel
            st.pop_back();
        st.push_back(pts[i]);
    }

    if (!include_collinear && st.size() == 2 && st[0] == st[1])
        st.pop_back();

    pts = st;
}

```

<div style="page-break-after: always;"></div>

### Interseccao

```cpp

// Requer: Template (Vetor, EPS, ccw)
bool onSegment(Vetor a, Vetor b, Vetor p) {//p pertence ao segmento ab
    return ccw(a,b,p)==0 &&
       min(a.x,b.x)-EPS <= p.x && p.x <= max(a.x,b.x)+EPS &&
       min(a.y,b.y)-EPS <= p.y && p.y <= max(a.y,b.y)+EPS;
}

bool intersect(Vetor a, Vetor b, Vetor c, Vetor d) {//segmentos ab e cd se intersectam
    int d1 = ccw(a,b,c), d2 = ccw(a,b,d);
    int d3 = ccw(c,d,a), d4 = ccw(c,d,b);
    if(d1*d2<0 && d3*d4<0) return true;
    if(d1==0 && onSegment(a,b,c)) return true;
    if(d2==0 && onSegment(a,b,d)) return true;
    if(d3==0 && onSegment(c,d,a)) return true;
    if(d4==0 && onSegment(c,d,b)) return true;
    return false;
}

```

<div style="page-break-after: always;"></div>

### Poligono

```cpp

// Requer: Template (Vetor, ccw)
double area(vector<Vetor>& poly){//Funciona para poligonos convexos ou nao
    double s=0;
    int n=poly.size();
    for(int i=0;i<n;i++)
        s += poly[i] ^ poly[(i+1)%n];
    return fabs(s)/2;
}

bool inPolygon(vector<Vetor> polygon,Vetor p){//Poligono convexo em anti-horario, O(log n)
    Vetor piv=polygon[0];
    if(p==piv) return true;
    if(ccw(p,piv,polygon[1])<0) return false;
    if(ccw(polygon[polygon.size()-1],piv,p)<0) return false;
    int l=1,r=polygon.size()-1;
    while(l+1<r){
        int mid = (l+r)>>1;
        if(ccw(piv,polygon[mid],p)>=0) l=mid;
        else r=mid;
    }
    return ccw(polygon[l],polygon[r],p)>=0;
}

```

<div style="page-break-after: always;"></div>

### Reta

```cpp

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

```

<div style="page-break-after: always;"></div>

### Template

```cpp

typedef long long int ll;
const double EPS = 1e-9;
const double PI = acos(-1);

// Vetor: +(add) -(sub) *(double)=escala *(Vetor)=dot ^=cross ~=modulo dist(p,q)=dist ao segmento pq
// cmpPolar: usar com sort + Vetor::pivot para ordenar por angulo polar
class Vetor{
    public:
    double x,y;//Pode ser trocado por long long int dependendo da questao
    Vetor operator+(Vetor q) const{
        return {this->x+q.x,this->y+q.y};
    }
    Vetor operator-(Vetor q) const{
        return {this->x-q.x,this->y-q.y};
    }
    Vetor operator*(double k) const{//Escalar * vetor
        return {this->x*k,this->y*k};
    }
    Vetor operator/(double k) const{//Vetor / escalar
        return {this->x/k,this->y/k};
    }
    double operator*(Vetor q) const{//Produto escalar
        return this->x*q.x+this->y*q.y;
    }
    double operator^(Vetor q) const{//Produto vetorial
        return this->x*q.y-q.x*this->y;
    }
    double operator~() const{//Modulo / distancia da origem
        return sqrt((*this)*(*this));
    }
    double dist(Vetor p,Vetor q) const{//Distancia de *this ao segmento pq
        Vetor r=*this;
        if((q-p)*(r-p)<=0||(p-q)*(r-q)<=0) return min(~(r-p),~(r-q));
        else return abs(((r-p)^(q-p))/(~(q-p)));
    }
    bool operator==(const Vetor p)const {
        return fabs(this->x - p.x) < EPS && fabs(this->y - p.y) < EPS;
    }
    bool operator!=(const Vetor p)const {
        return !(*this==p);
    }
    bool operator<(const Vetor p)const {
        if(this->x!=p.x) return this->x<p.x;
        return this->y<p.y;
    }
    static Vetor pivot;
    static bool cmpPolar(const Vetor& a, const Vetor& b){//Ordena por angulo polar a partir de pivot
        Vetor A = {a.x - pivot.x, a.y - pivot.y};
        Vetor B = {b.x - pivot.x, b.y - pivot.y};
        ll cross = A ^ B;
        if(cross == 0)
            return (A.x*A.x + A.y*A.y) < (B.x*B.x + B.y*B.y);
        return cross < 0; // <0 horario, >0 anti-horario
    }
};
Vetor Vetor::pivot = {0,0};

Vetor centroide(const vector<Vetor>& pts){
    ll sx=0, sy=0;
    for(auto &p: pts){ sx += p.x; sy += p.y; }
    return {sx/(ll)pts.size(), sy/(ll)pts.size()};
}

class Circle{
    public:
    Vetor o;
    double r;
    Circle(Vetor p,double raio): o(p),r(raio){ }
    Circle(double x,double y,double raio): o({x,y}),r(raio){ }
};

int ccw(Vetor a, Vetor b, Vetor c) {//1=esquerda, -1=direita, 0=colinear
    double cross = (b-a)^(c-a);
    if(fabs(cross) < EPS) return 0;
    return cross > 0 ? 1 : -1;
}

```

<div style="page-break-after: always;"></div>

## Grafos

### BFS

```cpp
//O(n+a) sendo a o numero de arestas
#define NMAX 1000
vector<int> adj[NMAX];
int vis[NMAX], anterior[NMAX];

int BFS(int ini, int fim, int tam){ // 1 se tiver caminho, 0 caso nao
    queue<int> fila;
    for (int i = 0; i < tam; i++) vis[i] = 0;
    fila.push(ini); vis[ini] = 1; anterior[ini] = -1;
    while (!fila.empty()){
        int u = fila.front(); fila.pop();
        for (int v : adj[u]){
            if (!vis[v]){
                if (v == fim){ anterior[v] = u; return 1; }
                fila.push(v); anterior[v] = u; vis[v] = 1;
            }
        }
    }
    return 0;
}

```

<div style="page-break-after: always;"></div>

### Bellman-Ford

```cpp
#define INF 0x3F3F3F3F
#define NMAX 100
// Caminho minimo com aresta negativa
// O(V*E)
struct Edge { int u, v, w; };
vector<Edge> edges;
int custo[NMAX], anterior[NMAX];

// Retorna true se existir ciclo negativo
bool bellmanFord(int s, int n){
    for (int i = 0; i < n; i++){
        custo[i] = INF;
        anterior[i] = -1;
    }
    custo[s] = 0;

    for (int k = 0; k < n - 1; k++)
        for (auto& e : edges)
            if (custo[e.u] != INF && custo[e.v] > custo[e.u] + e.w){
                custo[e.v] = custo[e.u] + e.w;
                anterior[e.v] = e.u;
            }

    for (auto& e : edges)
        if (custo[e.u] != INF && custo[e.v] > custo[e.u] + e.w)
            return true; // Ciclo negativo

    return false;
}

```

<div style="page-break-after: always;"></div>

### DFS

```cpp
// deteccao de ciclos, tempo de entrada e saida, nivel caso seja arvore

const int N=11;

int cnt=1;
int in[N],out[N],depth[N];//inicializa in e out com 0
int dis[N],low[N];
vector<int> adj[N];
set<int> artPoints;
set<pair<int,int>> bridges;


int DFS(int v,int nivel){//DFS com tempo de in e out, detecta ciclo e permite dizer se vertice faz parte da subarvore do outro
    in[v]=cnt++;
    depth[v]=nivel;
    for(int i=0;i<adj[v].size();i++){
        if(in[adj[v][i]]&&!out[adj[v][i]]) return 1;
        if(!in[adj[v][i]]){
            if(DFS(adj[v][i],nivel+1)) return 1;
        }
    }
    out[v]=cnt++;
    return 0;
}

void dfs(int v,int par){//dfs com lowlink e dis para deteccao de pontes e pontos de articulacao
    bool parentEdge=false;
    int children=0;
    dis[v] = low[v] = cnt++;
    for(auto i:adj[v]){
        if(i==par&&!parentEdge){
            parentEdge=true;
            continue;
        }
        if(!dis[i]){
            children++;
            dfs(i,v);
            low[v] = min(low[v], low[i]);
            if (par!=-1 && low[i] >= dis[v]) artPoints.insert(v);
            if (low[i] > dis[v]) bridges.insert({min(v,i),max(v,i)});
        }else low[v] = min(low[v], dis[i]);
    }
    if (par==-1 && children>1) artPoints.insert(v);
}
```

<div style="page-break-after: always;"></div>

### Dijkstra

```cpp
// Menor caminho
//O(e+nlogn) sendo e=arestas
#define SIZE 1000
#define INF 0x3f3f3f3f
typedef pair<int, int> ii;
typedef vector<ii> vii;
typedef vector<int> vi;

vii adj[SIZE];
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

### Dinic

```cpp
//O(EV^2)
struct FlowEdge {
    int v, u;
    long long cap, flow = 0;
    FlowEdge(int v, int u, long long cap) : v(v), u(u), cap(cap) {}
};

struct Dinic {
    const long long flow_inf = 1e18;
    vector<FlowEdge> edges;
    vector<vector<int>> adj;
    int n, m = 0;
    int s, t;
    vector<int> level, ptr;
    queue<int> q;

    Dinic(int n, int s, int t) : n(n), s(s), t(t) {
        adj.resize(n);
        level.resize(n);
        ptr.resize(n);
    }

    void add_edge(int v, int u, long long cap) {
        edges.emplace_back(v, u, cap);
        edges.emplace_back(u, v, 0);
        adj[v].push_back(m);
        adj[u].push_back(m + 1);
        m += 2;
    }

    bool bfs() {
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (int id : adj[v]) {
                if (edges[id].cap == edges[id].flow)
                    continue;
                if (level[edges[id].u] != -1)
                    continue;
                level[edges[id].u] = level[v] + 1;
                q.push(edges[id].u);
            }
        }
        return level[t] != -1;
    }

    long long dfs(int v, long long pushed) {
        if (pushed == 0)
            return 0;
        if (v == t)
            return pushed;
        for (int& cid = ptr[v]; cid < (int)adj[v].size(); cid++) {
            int id = adj[v][cid];
            int u = edges[id].u;
            if (level[v] + 1 != level[u])
                continue;
            long long tr = dfs(u, min(pushed, edges[id].cap - edges[id].flow));
            if (tr == 0)
                continue;
            edges[id].flow += tr;
            edges[id ^ 1].flow -= tr;
            return tr;
        }
        return 0;
    }

    long long flow() {
        long long f = 0;
        while (true) {
            fill(level.begin(), level.end(), -1);
            level[s] = 0;
            q.push(s);
            if (!bfs())
                break;
            fill(ptr.begin(), ptr.end(), 0);
            while (long long pushed = dfs(s, flow_inf)) {
                f += pushed;
            }
        }
        return f;
    }
};

int main(){
    int n,m,src,snk;
    cin>>n>>m>>src>>snk;
    Dinic ans(n,src,snk);
    while(m--){
        int a,b,c;
        cin>>a>>b>>c;
        ans.add_edge(a,b,c);
        //ans.add_edge(b,a,c) se bidirecional
    }
    cout<<ans.flow()<<endl;
    return 0;
}
```

<div style="page-break-after: always;"></div>

### Edmonds-Karp

```cpp
//O(Ef) E=num de arestas, f=fluxomax
struct Edge{
    int to;
    int cap;
};

const int N=112,INF=0x3f3f3f3f;

int vis[N],anterior[N],aresta[N];
vector<int> g[N];
vector<Edge> edges;
set<int> s;
//Funcao auxiliar para criar arestas
void addEdge(int a,int b,int c){
    g[a].push_back(edges.size());
    edges.push_back(Edge{b,c});
    g[b].push_back(edges.size());
    edges.push_back(Edge{a,c}); // Para arestas direcionadas use capacidade 0 na aresta reversa
}
//Roda BFS para pegar o caminho
int BFS(int ini, int fim){ // 1 se tiver caminho, 0 caso nao
    int i;
    queue<int> fila;
    memset(vis,0,sizeof(vis));
    fila.push(ini); vis[ini] = 1; anterior[ini] = -1;
    while (!fila.empty()){
        int a=fila.front();fila.pop();
        for(auto i:g[a]){
            int b=edges[i].to,c=edges[i].cap;
            if(!vis[b]&&c){
                anterior[b]=a;
                aresta[b]=i;
                if(b==fim) return 1;
                fila.push(b);
                vis[b]=1;
            }
        }
    }
    return 0;
}

int fluxoMaximo(int ini, int fim){
    int u, v;
    int fluxo = 0;
    int bot;
    while (BFS(ini, fim)){
        bot = INF;
        for (v = fim; v != ini; v = anterior[v]){
            u = aresta[v];
            bot = min(bot, edges[u].cap);
        }
        for (v = fim; v != ini; v = anterior[v]){
            u = aresta[v];
            edges[u].cap -= bot; edges[u^1].cap += bot; 
        }
        fluxo += bot;
    }
    return fluxo;
}
//Funcao que determina arestas que participam do minCut
int minCut(int tam){
    for(int i=0;i<tam;i++){
        if(vis[i]){
            for(auto j:g[i]){
                int b=edges[j].to;
                if(!vis[b]){
                    s.insert(j/2+1);
                }
            }
        }
    }
    return s.size();
}

```

<div style="page-break-after: always;"></div>

### Floyd-Warshall

```cpp
// Deteccao de ciclo negativo e caminho minimo para qualquer u, v
//O(n^3)
#define NMAX 1000
#define INF 0x3f3f3f3f
int m[NMAX][NMAX], custo[NMAX][NMAX];

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
            if (custo[i][i]< 0)
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
#define NMAX 1000
int pai[NMAX], rnk[NMAX];

int find(int u){
    return pai[u] = (pai[u] == u ? u : find(pai[u]));
}

void merge(int u, int v){
    u = find(u); v = find(v);
    if(rnk[u] > rnk[v]) pai[v] = u; else pai[u] = v;
    if(rnk[u] == rnk[v]) rnk[v]++;
}

// edges: vetor de (peso, u, v)
int kruskall(int n, vector<tuple<int,int,int>>& edges){
    for (int i = 0; i < n; i++){
        pai[i] = i; rnk[i] = 0;
    }
    sort(edges.begin(), edges.end());
    int res = 0;
    for (auto& [w, u, v] : edges)
        if (find(u) != find(v)){
            res += w; merge(u, v);
        }
    return res;
}

```

<div style="page-break-after: always;"></div>

### Least-Common-Ancestor

```cpp

const int N=212345;
const int LN=20;
int depth[N],memo[LN][N];
vector<int> filhos[N];
//O(nlogn)
void dfs(int s, int par){
    memo[0][s]=par;
    for(int i=0;i<filhos[s].size();i++){
        int f=filhos[s][i];
        if(f==par) continue;
        depth[f]=depth[s]+1;
        dfs(f,s);
    }
}
void compute(int n){
    for(int passos=1;passos<LN;passos++){
        for(int i=0;i<n;i++){
            memo[passos][i]=memo[passos-1][memo[passos-1][i]];
        }
    }
}
int walk(int s,int dis){
    int jump=0;
    while (dis){
        if(dis&1) s=memo[jump][s];
        jump++;
        dis>>=1;
    }
    return s;
}
int lca(int x,int y){
    if(depth[x]<depth[y]) swap(x,y);
    x=walk(x,depth[x]-depth[y]);
    if(x==y) return x;
    for(int k=LN-1;k>=0;k--){
        if(memo[k][x]!=memo[k][y]){
            x=memo[k][x];
            y=memo[k][y];
        }
    }
    return memo[0][x];
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
vector<int> comp[NMAX];
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
            comp[c].push_back(w); ins[w] = 0; s.pop();
        }
        w = s.top();
        comp[c].push_back(w); ins[w] = 0; s.pop();
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

### TopologicalSort

```cpp
//O(n+m)
const int N = 112345;
vector<int> adj[N];
bool vis[N],in[N];
vector<int> ans;

bool dfs(int v) {
    in[v]=vis[v]=true;
    for (int u : adj[v]) {
        if (!vis[u]) {
            if(dfs(u)) return true;
        }else if(in[u]) return true;
    }
    in[v]=false;
    ans.push_back(v);
    return false;
}

bool topological_sort(int n) {//Ordena topologicamente e caso seja impossivel retorna false
    ans.clear();
    for (int i = 0; i < n; ++i) {
        if (!vis[i]) {
            if(dfs(i)) return false;
        }
    }
    reverse(ans.begin(), ans.end());
    return true;
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

### Busca_Binaria

```cpp

const int INF=0x3f3f3f3f;

bool test(int n,int valor);//Funcao de teste no vetor

int bb(int n){//Busca Binaria que procura um valor que compra os requisitos da funcao teste
    int lmin=0,lmax=INF;
    int res=-1;
    while(lmin<=lmax){
        int mid=(lmin+lmax)/2;
        if(test(n,mid)){
            res=mid;
            lmin=mid+1;
        }
        else lmax=mid-1;
    }
    return res;
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
    vector<vector<double>> mat(n, vector<double>(n));
    double res=1;
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
    vector<double> aux(n);
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
const int N=11234567;
int fprimos[N]; // inicializa com 0
//p=primo !fprimos[p]||fprimos[p]==p
//p!=primo fprimos[p]&&fprimos[p]!=p
void crivo(){
    fprimos[0] = fprimos[1] = 1;
    for(int i = 2; i*i < N; i++){
        if(fprimos[i] == 0) fprimos[i] = i;
        for (int j = i*i; j < N; j+=i) fprimos[j] = fprimos[i];
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

### Inverso-Modular

```cpp
//Inverso Modular para divisão com MOD
//O(logMOD) precisa de gcd(n,MOD)=1 MOD precisa ser primo
const long long MOD=998244353;

long long mod_inv(long long num){return fpow(num,MOD-2);}

long long fpow(long long base, long long power) {
    long long result = 1;
    while(power > 0) {
        if(power&1) result = (result*base)%MOD ;
        base = (base * base)%MOD;
        power>>=1;
    }
    return result;
}

```

<div style="page-break-after: always;"></div>

### Segmented_Sieve

```cpp

// computa os primos num intervalo L..R
// O((R-L+1)loglogR + sqrt(R)loglog(sqrt(R)))
vector<char> segmentedSieve(long long L, long long R)
{
    // generate all primes up to sqrt(R)
    long long lim = sqrt(R);
    vector<char> mark(lim + 1, false);
    vector<long long> primes;
    for (long long i = 2; i <= lim; ++i)
    {
        if (!mark[i])
        {
            primes.emplace_back(i);
            for (long long j = i * i; j <= lim; j += i)
                mark[j] = true;
        }
    }

    vector<char> isPrime(R - L + 1, true);
    for (long long i : primes)
        for (long long j = max(i * i, (L + i - 1) / i * i); j <= R; j += i) // i*i decreases recalculations
            isPrime[j - L] = false;
    if (L == 1)
        isPrime[0] = false;
    return isPrime;
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
#define NMAX 1000
int BIT[NMAX];

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
    for (i = 0; i < tam; i++) updateBIT(tam, i, vet[i]);
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
#define NMAX 1000
int peso[NMAX], valor[NMAX];

// com repeticao
int knapSackRep(int W, int n){
    int memo[W+1];
    memset(memo, 0, sizeof(memo));
    for (int w = 0; w <= W; w++){
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

// versao com matriz. Permite reconstruir os itens escolhidos.
int mat[NMAX][NMAX];

int knapSackMat(int W, int n) {
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

### Persistent-Segment-Tree

```cpp

#define N 100

struct node
{
    int val;
    node* left, *right;

    node() {}
    node(node* l, node* r, int v)
    {
        left = l;
        right = r;
        val = v;
    }
};

int arr[N];

node* version[N];

//O(nlogn)
void build(node* n,int low,int high)
{
    if (low==high)
    {
        n->val = arr[low];
        return;
    }
    int mid = (low+high) / 2;
    n->left = new node(NULL, NULL, 0);
    n->right = new node(NULL, NULL, 0);
    build(n->left, low, mid);
    build(n->right, mid+1, high);
    n->val = n->left->val + n->right->val;
}

/*
 * Time Complexity : O(logn)
 * Space Complexity : O(logn)  */
void upgrade(node* prev, node* cur, int low, int high, int idx, int value)
{
    if (idx > high or idx < low or low > high)
        return;

    if (low == high)
    {
        // modification in new version
        cur->val = value;
        return;
    }
    int mid = (low+high) / 2;
    if (idx <= mid)
    {
        // link to right child of previous version
        cur->right = prev->right;

        // create new node in current version
        cur->left = new node(NULL, NULL, 0);

        upgrade(prev->left,cur->left, low, mid, idx, value);
    }
    else
    {
        // link to left child of previous version
        cur->left = prev->left;

        // create new node for current version
        cur->right = new node(NULL, NULL, 0);

        upgrade(prev->right, cur->right, mid+1, high, idx, value);
    }

    // calculating data for current version
    cur->val = cur->left->val + cur->right->val;
}

int query(node* n, int low, int high, int l, int r)
{
    if (l > high or r < low or low > high)
       return 0;
    if (l <= low and high <= r)
       return n->val;
    int mid = (low+high) / 2;
    int p1 = query(n->left,low,mid,l,r);
    int p2 = query(n->right,mid+1,high,l,r);
    return p1+p2;
}

int main(int argc, char const *argv[])
{
    int A[] = {1,2,3,4,5};
    int n = sizeof(A)/sizeof(int);

    for (int i=0; i<n; i++) 
       arr[i] = A[i];

    // creating Version-0
    node* root = new node(NULL, NULL, 0);
    build(root, 0, n-1);

    // storing root node for version-0
    version[0] = root;

    // upgrading to version-1
    version[1] = new node(NULL, NULL, 0);
    upgrade(version[0], version[1], 0, n-1, 4, 1);

    // upgrading to version-2
    version[2] = new node(NULL, NULL, 0);
    upgrade(version[1],version[2], 0, n-1, 2, 10);
    
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
const int N=1123;
int vet[N];
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

### SparseTable

```cpp
//O(nlogn) construcao O(1)ouO(logn) para querys
const int K=21;
const int N=112345;
int st[K + 1][N];
int lg[N + 1];
vector<int> v;

int f(int a,int b){//Funcao das querys
    return min(a,b);
}
void build(){
    copy(v.begin(), v.end(), st[0]);
    lg[1] = 0;
    for (int i = 2; i <= N; i++)
        lg[i] = lg[i/2] + 1;

    for (int i = 1; i <= K; i++)
        for (int j = 0; j + (1 << i) <= N; j++)
            st[i][j] = f(st[i - 1][j], st[i - 1][j + (1 << (i - 1))]);
}

int queryI(int L,int R){//Funcao Idempotente: O(1)
    int i = lg[R - L + 1];
    return f(st[i][L], st[i][R - (1 << i) + 1]);
}

long long query(int L,int R){//Exemplo para soma: O(logn)
    long long sum = 0;
    for (int i = K; i >= 0; i--) {
        if ((1 << i) <= R - L + 1) {
            sum += st[i][L];
            L += 1 << i;
        }
    }
    return sum;
}

```

<div style="page-break-after: always;"></div>

### SubsetSum

```cpp
//O(n*sum) nao funciona com numero negativo
vector<int> nums;

bool subsetsum(int sum) {
    vector<bool> prev(sum+1,false),cur(sum+1,false);
    prev[0]=true;
    for(auto i:nums){
        for(int j=0;j<=sum;j++){
            if(j<i) cur[j]=prev[j];
            else cur[j]=(prev[j]||prev[j-i]);
        }
        prev=cur;
    }
    return prev[sum];
}

```

<div style="page-break-after: always;"></div>

### Sum-Array-2D

```cpp
//Cria uma matriz de soma, para verificar a soma num quadrado usar a funcao soma, pode ser modificado para retangulo
//O(n^2)
const int N=1123;
int prefix[N][N],vet[N][N];
void monta(int n,int m){
  memset(prefix,0,sizeof(prefix));
  for(int i=1;i<=n;i++){
    for(int j=1;j<=m;j++) prefix[i][j]=prefix[i-1][j]+prefix[i][j-1]-prefix[i-1][j-1]+vet[i-1][j-1];
  }
}

int soma(int x,int y,int lado){
  return prefix[x+lado][y+lado]-prefix[x+lado][y-1]-prefix[x-1][y+lado]+prefix[x-1][y-1];
}

```

<div style="page-break-after: always;"></div>

## String

### AhoCorasick

```cpp
//O(mk) m=tamanho total das strings, k=tamanho do alfabeto
// Chamar build() depois de adicionar todas as strings
const int K = 26;

struct Vertex {
    int next[K];
    int go[K];
    bool output = false;
    int link = 0;

    Vertex() {
        fill(begin(next), end(next), -1);
        fill(begin(go), end(go), 0);
    }
};

vector<Vertex> t(1);

void add_string(string const& s) {
    int v = 0;
    for (char ch : s) {
        int c = ch - 'a';
        if (t[v].next[c] == -1) {
            t[v].next[c] = t.size();
            t.emplace_back();
        }
        v = t[v].next[c];
    }
    t[v].output = true;
}

void build() {
    queue<int> q;
    for (int c = 0; c < K; c++) {
        if (t[0].next[c] == -1) {
            t[0].go[c] = 0;
        } else {
            t[0].go[c] = t[0].next[c];
            q.push(t[0].next[c]);
        }
    }
    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int c = 0; c < K; c++) {
            if (t[v].next[c] != -1) {
                int u = t[v].next[c];
                t[u].link = t[t[v].link].go[c];
                t[v].go[c] = u;
                q.push(u);
            } else {
                t[v].go[c] = t[t[v].link].go[c];
            }
        }
    }
}

```

<div style="page-break-after: always;"></div>

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

### Manacher

```cpp
//Funcao para Palindromos
const int S=200010;//Tamanho maximo da palavra

int odd[S];//Guarda o maior tamanho de palindromo centrado na posicao i, tamanho = 2*odd[i]-1
int even[S];//Guarda o maior tamanho de palindromo centrado nas posicoes i e i-1, tamanho = 2*even[i]

void manacher(string s) {//O(n)
    string t;
    for(auto c: s) {
        t += string("#") + c;
    }
    s = t + "#";
    int n = s.size();
    s = "$" + s + "^";
    vector<int> p(n + 2);
    int l = 0, r = 1;
    for(int i = 1; i <= n; i++) {
        p[i] = min(r - i, p[l + (r - i)]);
        while(s[i - p[i]] == s[i + p[i]]) {
            p[i]++;
        }
        if(i + p[i] > r) {
            l = i - p[i], r = i + p[i];
        }
    }
    for(int i=0;i<n;i++){
        if(i&1) odd[(i-1)/2]=p[i+1]/2;
        else even[i/2]=(p[i+1]-1)/2;
    }
}

```

<div style="page-break-after: always;"></div>

### Prefix-Function

```cpp

const int N=1123456;
int pi[N];//Caso for acessar tomar cuidado pois o valor dentro eh o tamanho
int ans[N];//Prefixo de tamanho i aparece ans[i] vezes
//Usado para contar prefixo, procurar substring e contar substrings
void prefix_function(string s){//Guarda o tamanho da string que eh sufixo e prefixo simultaneamente
    //A letra i eh o final de um sufixo e prefixo de tamanho pi[i]
    pi[0]=0;
    for(int i=1;i<s.length();i++){
        int j=pi[i-1];
        while(j>0&&s[j]!=s[i]) j=pi[j-1];
        if(s[j]==s[i]) pi[i]=j+1;
        else pi[i]=0; 
    }
}

void conta_prefixos(string s){//Conta quantas vezes cada prefixo aparece na string
    int n=s.length();
    prefix_function(s);
    for (int i = 0; i < n; i++)
        ans[pi[i]]++;
    for (int i = n-1; i > 0; i--)
        ans[pi[i-1]] += ans[i];
    for (int i = 0; i <= n; i++)
        ans[i]++;
}

```

<div style="page-break-after: always;"></div>

### StringHashing

```cpp

const int p = 31;
const int m = 1000000009;

long long compute_hash(string const& s) {//Calcula hash em O(S) mas permite comparações entre strings em O(1)
    long long hash_value = 0;
    long long p_pow = 1;
    for (char c : s) {
        hash_value = (hash_value + (c - 'a' + 1) * p_pow) % m;
        p_pow = (p_pow * p) % m;
    }
    return hash_value;
}

vector<vector<int>> group_identical_strings(vector<string> const& s) {//Utiliza hash para agrupar strings iguais
    int n = s.size();
    vector<pair<long long, int>> hashes(n);
    for (int i = 0; i < n; i++)
        hashes[i] = {compute_hash(s[i]), i};

    sort(hashes.begin(), hashes.end());

    vector<vector<int>> groups;
    for (int i = 0; i < n; i++) {
        if (i == 0 || hashes[i].first != hashes[i-1].first)
            groups.emplace_back();
        groups.back().push_back(hashes[i].second);
    }
    return groups;
}

vector<int> rabin_karp(string const& s, string const& t) {//Conta quantas ocorrencias de s existe no texto t O(T+S)
    int S = s.size(), T = t.size();

    vector<long long> p_pow(max(S, T)); 
    p_pow[0] = 1; 
    for (int i = 1; i < (int)p_pow.size(); i++) 
        p_pow[i] = (p_pow[i-1] * p) % m;

    vector<long long> h(T + 1, 0); 
    for (int i = 0; i < T; i++)
        h[i+1] = (h[i] + (t[i] - 'a' + 1) * p_pow[i]) % m; 
    long long h_s = 0; 
    for (int i = 0; i < S; i++) 
        h_s = (h_s + (s[i] - 'a' + 1) * p_pow[i]) % m; 

    vector<int> occurrences;
    for (int i = 0; i + S - 1 < T; i++) {
        long long cur_h = (h[i+S] + m - h[i]) % m;
        if (cur_h == h_s * p_pow[i] % m)
            occurrences.push_back(i);
    }
    return occurrences;
}
```

<div style="page-break-after: always;"></div>

### SuffixAutomaton

```cpp
//O(n)
struct state {
    int len, link;
    map<char, int> next;
};

const int N = 112345;
state st[N * 2];
int sz, last;

void sa_init() {
    st[0].len = 0;
    st[0].link = -1;
    st[0].next.clear();
    sz = 1;
    last = 0;
}

void sa_extend(char c) {
    int cur = sz++;
    st[cur].len = st[last].len + 1;
    int p = last;
    while (p != -1 && !st[p].next.count(c)) {
        st[p].next[c] = cur;
        p = st[p].link;
    }
    if (p == -1) {
        st[cur].link = 0;
    } else {
        int q = st[p].next[c];
        if (st[p].len + 1 == st[q].len) {
            st[cur].link = q;
        } else {
            int clone = sz++;
            st[clone].len = st[p].len + 1;
            st[clone].next = st[q].next;
            st[clone].link = st[q].link;
            while (p != -1 && st[p].next[c] == q) {
                st[p].next[c] = clone;
                p = st[p].link;
            }
            st[q].link = st[cur].link = clone;
        }
    }
    last = cur;
}

long long get_diff_strings(){//Total de substrings diferentes
    long long tot = 0;
    for(int i = 1; i < sz; i++) {
        tot += st[i].len - st[st[i].link].len;
    }
    return tot;
}

long long get_tot_len_diff_substrings() {//Soma dos tamanhos de substrings diferentes
    long long tot = 0;
    for(int i = 1; i < sz; i++) {
        long long shortest = st[st[i].link].len + 1;
        long long longest = st[i].len;

        long long num_strings = longest - shortest + 1;
        long long cur = num_strings * (longest + shortest) / 2;
        tot += cur;
    }
    return tot;
}

string lcs (string S, string T) {//Longest Common Substring
    sa_init();
    for (int i = 0; i < S.size(); i++)
        sa_extend(S[i]);

    int v = 0, l = 0, best = 0, bestpos = 0;
    for (int i = 0; i < T.size(); i++) {
        while (v && !st[v].next.count(T[i])) {
            v = st[v].link ;
            l = st[v].len;
        }
        if (st[v].next.count(T[i])) {
            v = st [v].next[T[i]];
            l++;
        }
        if (l > best) {
            best = l;
            bestpos = i;
        }
    }
    return T.substr(bestpos - best + 1, best);
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
