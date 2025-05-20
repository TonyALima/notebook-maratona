#include <bits/stdc++.h>
using namespace std;
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

void convex_hull(vector<pair<pt,int>>& a, bool include_collinear = false) {
    pair<pt,int> p0 = *min_element(a.begin(), a.end(), [](auto a, auto b) {
        return make_pair(a.first.y, a.first.x) < make_pair(b.first.y, b.first.x);
    });
    sort(a.begin(), a.end(), [&p0](const auto& a, const auto& b) {
        int o = orientation(p0.first, a.first, b.first);
        if (o == 0)
            return (p0.first.x-a.first.x)*(p0.first.x-a.first.x) + (p0.first.y-a.first.y)*(p0.first.y-a.first.y)
                < (p0.first.x-b.first.x)*(p0.first.x-b.first.x) + (p0.first.y-b.first.y)*(p0.first.y-b.first.y);
        return o < 0;
    });
    if (include_collinear) {
        int i = (int)a.size()-1;
        while (i >= 0 && collinear(p0.first, a[i].first, a.back().first)) i--;
        reverse(a.begin()+i+1, a.end());
    }

    vector<pair<pt,int>> st;
    for (int i = 0; i < (int)a.size(); i++) {
        while (st.size() > 1 && !cw(st[st.size()-2].first, st.back().first, a[i].first, include_collinear))
            st.pop_back();
        st.push_back(a[i]);
    }

    if (include_collinear == false && st.size() == 2 && st[0] == st[1])
        st.pop_back();

    a = st;
}

bool cmp(pair<pt,int> x,pair<pt,int> y){
    return x.second<y.second;
}

int main(){
    int n;
    vector<pair<pt,int>> v;
    cin>>n;
    for(int i=1;i<=n;i++){
        double x,y;
        cin>>x>>y;
        pt a;
        a.x=x;a.y=y;
        v.push_back(make_pair(a,i));
    }
    convex_hull(v,true);
    sort(v.begin(),v.end(),cmp);
    cout<<v[0].second;
    for(int i=1;i<v.size();i++) cout<<" "<<v[i].second;
    cout<<endl;
}
