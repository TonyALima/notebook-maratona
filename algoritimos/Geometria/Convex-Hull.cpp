#include <bits/stdc++.h>
using namespace std;
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
