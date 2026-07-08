#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <queue>
#include <tuple>
#include <stack>

using namespace std;

const double EPS=1e-9;

struct half_ed;
struct face;

struct point {
    double x, y;
    half_ed* rand_one_ed;  
    bool is_hole;          
    bool is_rightmost;     
    int idx;               

    point(double x_i=0, double y_i=0)
        : x(x_i), y(y_i),rand_one_ed(nullptr),
          is_hole(false),is_rightmost(false),idx(-1) {}
};

struct half_ed {
    point* st;            
    half_ed* next;          
    half_ed* prev;          
    half_ed* reverse_ed;    
    face* incident_face; 

    half_ed()
        : st(nullptr),next(nullptr),prev(nullptr),
          reverse_ed(nullptr),incident_face(nullptr) {}
};

double getSignedArea(const vector<pair<double,double>>& poly) {
    double area=0;
    int n=(int)poly.size();
    for (int i=0;i<n;i++) {
        int j=(i+1) % n;
        area +=poly[i].first*poly[j].second
              - poly[j].first*poly[i].second;
    }
    return area / 2.0;
}

double cross(point* a, point* b, point* c) {
    return (b->x-a->x)*(c->y-a->y)-(b->y-a->y)*(c->x-a->x);
}

bool above(point* a, point* b) {
    if (fabs(a->y-b->y)>EPS) return a->y>b->y;
    return a->x<b->x;
}

enum Type { START, END, SPLIT, MERGE, REGULAR };

struct DCEL {
    vector<point*> points;
    vector<half_ed*> half_eds;
    vector<pair<int,int>> diagonals;
    
    vector<pair<int,int>> original_edges;
    vector<int> orig_prev;     
    vector<int> orig_next;     
    
    vector<tuple<int,int,int>> final_triangles;

    int n_original=0;  
    int h_count =0;  

    point* join_pt(double x, double y) {
        point* p=new point(x, y);
        p->idx=(int)points.size();
        points.push_back(p);
        return p;
    }

    void add_ring(const vector<point*>& ring_pts) {
        int n=(int)ring_pts.size();
        vector<half_ed*> edges(n);
        for (int i=0;i<n;i++) {
            edges[i]=new half_ed();
            edges[i]->st=ring_pts[i];
            ring_pts[i]->rand_one_ed=edges[i];
            half_eds.push_back(edges[i]);
        }
        for (int i=0;i<n;i++) {
            edges[i]->next=edges[(i+1)%n];
            edges[i]->prev=edges[(i-1+n)%n];
            original_edges.push_back({ring_pts[i]->idx, ring_pts[(i+1)%n]->idx});
        }
    }

    void add_diagonal(int u_idx, int v_idx) {
        diagonals.push_back({u_idx, v_idx});

        point* u=points[u_idx];
        point* v=points[v_idx];
        
        half_ed* he_u=u->rand_one_ed;
        half_ed* he_v=v->rand_one_ed;

        half_ed* d_uv=new half_ed();
        half_ed* d_vu=new half_ed();

        d_uv->st=u;
        d_vu->st=v;
        d_uv->reverse_ed=d_vu;
        d_vu->reverse_ed=d_uv;

        half_eds.push_back(d_uv);
        half_eds.push_back(d_vu);

        half_ed* he_u_prev=he_u->prev;
        half_ed* he_v_prev=he_v->prev;

        he_u_prev->next=d_uv; 
        d_uv->prev=he_u_prev;
        d_uv->next=he_v;      
        he_v->prev=d_uv;

        he_v_prev->next=d_vu; 
        d_vu->prev=he_v_prev;
        d_vu->next=he_u;      
        he_u->prev=d_vu;
    }

    void clear() {
        for (auto p : points)    delete p;
        for (auto e : half_eds)  delete e;
        points.clear();
        half_eds.clear();
        diagonals.clear();
        original_edges.clear();
        orig_prev.clear();
        orig_next.clear();
        final_triangles.clear();
        n_original=0;
        h_count =0;
    }

    vector<int> three_color() {
        int n=(int)points.size();
        vector<int> color(n,-1);
        if (final_triangles.empty()) return color;
        int T=(int)final_triangles.size();
        map<pair<int,int>, int> edge_to_tri;
        for (int i=0;i<T;i++) {
            auto [a, b, c]=final_triangles[i];
            edge_to_tri[{a,b}]=i;
            edge_to_tri[{b,c}]=i;
            edge_to_tri[{c,a}]=i;
        }

        vector<bool> visited(T, false);
        queue<int> bfs;
        auto propagate_triangle=[&](int ti) {
            auto [a, b, c]=final_triangles[ti];
            int verts[3]={a, b, c};
            for (int v : verts) {
                if (color[v]!=-1) continue;
                bool used[3]={false, false, false};
                for (int w : verts) {
                    if (w!=v && color[w]!=-1) used[color[w]]=true;
                }
                for (int col=0;col<3;col++) {
                    if (!used[col]) {color[v]=col;break;}
                }
            }
        };
        auto bfs_component=[&](int seed) {
            visited[seed]=true;
            bfs.push(seed);
            while(!bfs.empty()) {
                int ti=bfs.front();bfs.pop();
                auto [a, b, c]=final_triangles[ti];
                pair<int,int> edges[3]={{a,b},{b,c},{c,a}};
                for (auto& e : edges) {
                    auto it=edge_to_tri.find({e.second, e.first});
                    if (it==edge_to_tri.end()) continue;
                    int ni=it->second;
                    if (visited[ni]) continue;
                    visited[ni]=true;
                    propagate_triangle(ni);
                    bfs.push(ni);
                }
            }
        };

        {
            auto [a0, b0, c0]=final_triangles[0];
            color[a0]=0;color[b0]=1;color[c0]=2;
        }
        bfs_component(0);

        for (int i=0;i<T;i++) {
            if (visited[i]) continue;

            auto [a, b, c]=final_triangles[i];
            int verts[3]={a, b, c};
            bool used[3]={false, false, false};
            for (int v:verts) if (color[v]!=-1) used[color[v]]=true;
            for (int v:verts) {
                if (color[v]!=-1) continue;
                for (int col=0;col<3;col++) {
                    if (!used[col]) { color[v]=col;used[col]=true;break;}
                }
            }
            bfs_component(i);
        }

        return color;
    }
};

DCEL* g_dcel=nullptr;
vector<point*>* g_points=nullptr;
double sweepline=0.0;

#define points (*g_points)
#define PUSH_DIAG(a, b) g_dcel->add_diagonal((a),(b))

double getX(int i) {
    point* a=points[i];
    point* b=points[g_dcel->orig_next[i]];
    if (fabs(a->y - b->y)<EPS) return min(a->x, b->x);
    return a->x+(b->x-a->x)*(sweepline-a->y)/(b->y-a->y);
}

struct cmp {
    using is_transparent=void;
    bool operator()(int a, int b) const {
        double x1=getX(a), x2=getX(b);
        if (fabs(x1-x2)>EPS) return x1<x2;
        return a<b;
    }
    bool operator()(int a, double x) const { return getX(a)<x-EPS;}
    bool operator()(double x, int b) const { return x<getX(b)-EPS;}
};

Type getType(int i) {
    point* cur=points[i];
    point* prev_pt=points[g_dcel->orig_prev[i]];
    point* next_pt=points[g_dcel->orig_next[i]];

    bool prevAbove=above(prev_pt, cur);
    bool nextAbove=above(next_pt, cur);
    double cr=cross(prev_pt, cur, next_pt);

    if (!prevAbove && !nextAbove) return (cr<-EPS) ? SPLIT:START;
    if ( prevAbove && nextAbove) return (cr<-EPS) ? MERGE:END;
    return REGULAR;
}

void make_monotone(int n_original) {
    int n=(int)points.size();

    vector<int> ord(n);
    iota(ord.begin(),ord.end(),0);

    sort(ord.begin(),ord.end(),[](int a, int b) {
        return above(points[a],points[b]);
    });

    vector<Type> v_type(n, REGULAR);
    for (int i=0;i<n_original;i++)
        v_type[i]=getType(i);

    set<int, cmp> T;
    unordered_map<int,int> helper;

    vector<set<int, cmp>::iterator> active_iters(n_original,T.end());
    vector<bool> in_T(n_original,false);

    for (int v:ord) {
        if(v>=n_original) continue;

        sweepline=points[v]->y;
        Type t=v_type[v];

        int pv=g_dcel->orig_prev[v];
        point* curPt =points[v];
        point* nextPt=points[g_dcel->orig_next[v]];

        if(t==START){
            auto it=T.insert(v).first;
            active_iters[v]=it;
            in_T[v]=true;
            helper[v]=v;
        }

        else if(t==END){
            if (helper.count(pv) && v_type[helper[pv]]==MERGE)
                PUSH_DIAG(v, helper[pv]);

            if (in_T[pv]) {
                T.erase(active_iters[pv]);
                in_T[pv]=false;
            }
            helper.erase(pv);
        }

        else if(t==SPLIT){
            auto it=T.lower_bound(v);
            if (it !=T.begin()) {
                --it;
                int e=*it;

                if (helper.count(e))
                    PUSH_DIAG(v, helper[e]);

                helper[e]=v;
            }

            auto it2=T.insert(v).first;
            active_iters[v]=it2;
            in_T[v]=true;
            helper[v]=v;
        }

        else if(t==MERGE){
            if (helper.count(pv) && v_type[helper[pv]]==MERGE)
                PUSH_DIAG(v, helper[pv]);
            if (in_T[pv]) {
                T.erase(active_iters[pv]);
                in_T[pv]=false;
            }
            helper.erase(pv);

            auto it=T.lower_bound(v);
            if (it !=T.begin()) {
                --it;
                int e=*it;
                if (helper.count(e) && v_type[helper[e]]==MERGE)
                    PUSH_DIAG(v, helper[e]);

                helper[e]=v;
            }
        }
        else{ 
            bool is_left=above(curPt, nextPt);
            if(is_left){
                if (helper.count(pv) && v_type[helper[pv]]==MERGE)
                    PUSH_DIAG(v, helper[pv]);
                if (in_T[pv]) {
                    T.erase(active_iters[pv]);
                    in_T[pv]=false;
                }
                helper.erase(pv);
                auto it=T.insert(v).first;
                active_iters[v]=it;
                in_T[v]=true;
                helper[v]=v;
            }
            else {
                auto it=T.lower_bound(v);
                if (it !=T.begin()) {
                    --it;
                    int e=*it;

                    if (helper.count(e) && v_type[helper[e]]==MERGE)
                        PUSH_DIAG(v, helper[e]);

                    helper[e]=v;
                }
            }
        }
    }
}
#undef points
#undef PUSH_DIAG

vector<vector<int>> get_monotone_polygon(DCEL& dcel) {
    int n=dcel.n_original;
    vector<vector<int>> adj(n);
    for(auto& edge : dcel.original_edges) {
        adj[edge.first].push_back(edge.second);
        adj[edge.second].push_back(edge.first);
    }
    set<pair<int, int>> unique_diags;
    for(auto& diag : dcel.diagonals) {
        int u=min(diag.first, diag.second);
        int v=max(diag.first, diag.second);
        if (u<n && v<n) unique_diags.insert({u, v});
    }
    for(auto& diag : unique_diags) {
        adj[diag.first].push_back(diag.second);
        adj[diag.second].push_back(diag.first);
    }
for (int i=0;i<n;i++) {
    point* p=dcel.points[i];

    sort(adj[i].begin(), adj[i].end());
    adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
    sort(adj[i].begin(), adj[i].end(), [&](int a, int b) {
        point* A=dcel.points[a];
        point* B=dcel.points[b];
        double cross=(A->x - p->x)*(B->y - p->y) -
                       (A->y - p->y)*(B->x - p->x);
        if (cross==0) {
            double da=(A->x - p->x)*(A->x - p->x) +
                        (A->y - p->y)*(A->y - p->y);
            double db=(B->x - p->x)*(B->x - p->x) +
                        (B->y - p->y)*(B->y - p->y);
            return da<db;
        }
        return cross>0;
    });
        adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
    }
    vector<unordered_map<int,int>> rank(n);
    for(int i=0;i<n;i++)
        for(int j=0;j<(int)adj[i].size();j++)
            rank[i][adj[i][j]]=j;
    map<pair<int,int>, bool> visited;
    vector<vector<int>> polys;
    for(int i=0;i<n;i++) {
        for(int nxt : adj[i]) {
            if(visited[{i, nxt}]) continue;
            vector<int> loop;
            int u=i, v=nxt;
            while(!visited[{u, v}]) {
                visited[{u, v}]=true;
                loop.push_back(u);
                int idx=rank[v][u];
                int next_idx=(idx - 1 + (int)adj[v].size()) % (int)adj[v].size();
                int new_v=adj[v][next_idx];
                u=v;v=new_v;
            }
            if(loop.size() >=3) {
                vector<pair<double,double>> pts;
                for(int idx : loop) pts.push_back({dcel.points[idx]->x, dcel.points[idx]->y});
                if(getSignedArea(pts)>EPS) polys.push_back(loop);
            }
        }
    }

    vector<vector<int>> valid_polys;
    for(auto& loop : polys) {
        bool has_outer_vertex=false;
        for(int idx : loop) {
            if(!dcel.points[idx]->is_hole) { has_outer_vertex=true;break;}
        }
        if(has_outer_vertex) {
            valid_polys.push_back(loop);
            continue;
        }
        set<int> ring_ids;
        for(int idx : loop) {
            int cur=idx, min_in_ring=idx;
            do {
                min_in_ring=min(min_in_ring, cur);
                cur=dcel.orig_next[cur];
            } while(cur !=idx && cur !=-1);
            ring_ids.insert(min_in_ring);
        }
        if(ring_ids.size()>1) {
            valid_polys.push_back(loop);
        }
    }
    return valid_polys;
}
static double orient_tri(point* a, point* b, point* c) {
    return (b->x - a->x) * (c->y - a->y) - (b->y - a->y) * (c->x - a->x);
}
static vector<int> get_chain(vector<point*>& poly) {
    int n=(int)poly.size();
    int top=0, bot=0;
    for (int i=1;i<n;i++) {
        if (poly[i]->y>poly[top]->y || (poly[i]->y==poly[top]->y && poly[i]->x<poly[top]->x)) top=i;
        if (poly[i]->y<poly[bot]->y || (poly[i]->y==poly[bot]->y && poly[i]->x>poly[bot]->x)) bot=i;
    }
    vector<int> chain(n, -1);
    vector<int> path1, path2;
    int i=top;
    while(true) { path1.push_back(i);if (i==bot) break;i=(i+1)%n;}
    i=top;
    while(true) { path2.push_back(i);if (i==bot) break;i=(i-1+n)%n;}
    double o=orient_tri(poly[top],poly[path1[1]],poly[path2[1]]);
    if (o>0) {
        for (int v:path1) if (v!=top && v!=bot) chain[v]=0;
        for (int v:path2) if (v!=top && v!=bot) chain[v]=1;
    } else {
        for (int v : path1) if (v !=top && v !=bot) chain[v]=1;
        for (int v : path2) if (v !=top && v !=bot) chain[v]=0;
    }
    chain[top]=0;chain[bot] =1;
    return chain;
}

void triangulate_monotone(DCEL& dcel, const vector<int>& face_indices) {
    int n=(int)face_indices.size();
    if (n<3) return;
    if (n==3) {
        dcel.final_triangles.push_back({face_indices[0], face_indices[1], face_indices[2]});
        return;
    }

    vector<point*> poly(n);
    vector<int>    origIdx(n);
    for (int i=0;i<n;i++) {
        poly[i] =dcel.points[face_indices[i]];
        origIdx[i]=face_indices[i];
    }

    vector<int> chain=get_chain(poly);
    vector<int> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) {
        if (poly[a]->y !=poly[b]->y) return poly[a]->y>poly[b]->y;
        if (chain[a] !=chain[b])     return chain[a]<chain[b];
        if (chain[a]==0)            return poly[a]->x>poly[b]->x;
        return poly[a]->x<poly[b]->x;
    });

    int first=order[0], last=order[n - 1];
    stack<int> S;
    S.push(order[0]);S.push(order[1]);
    int stack_chain=chain[order[1]];

    auto emit_triangle=[&](int a, int b, int c) {
        dcel.final_triangles.push_back({origIdx[a], origIdx[b], origIdx[c]});
    };

    for (int j=2;j<n - 1;j++) {
        int curr=order[j];
        if (chain[curr] !=stack_chain) {
            int prev_v=-1;
            vector<int> popped;
            while(S.size()>1) { popped.push_back(S.top());S.pop();}
            int bottom_v=S.top();

            for (int k=0;k<(int)popped.size();k++) {
                int v=popped[k];
                int u=(k + 1<(int)popped.size()) ? popped[k + 1] : bottom_v;
                emit_triangle(curr, u, v);
            }
            while(!S.empty()) S.pop();
            S.push(popped.back());S.push(curr);
            stack_chain=chain[curr];
        } else {
            int v=S.top();S.pop();
            while(!S.empty()) {
                int u=S.top();
                double o=orient_tri(poly[u], poly[v], poly[curr]);
                bool valid=(chain[curr]==0 && o>0) || (chain[curr]==1 && o<0);
                if (valid) {
                    emit_triangle(curr, u, v);
                    v=S.top();S.pop();
                } else break;
            }
            S.push(v);S.push(curr);
        }
    }

    int curr=last;
    int prev_v=-1;
    while(!S.empty()) {
        int v=S.top();S.pop();
        if (v==first) { prev_v=v;continue;}
        if (prev_v !=-1 && prev_v !=first) emit_triangle(curr, prev_v, v);
        prev_v=v;
    }
}

int main() {
    ifstream fin("input.txt");
    if (!fin.is_open()) {
        cerr << "Error: Could not open input.txt\n";
        return 1;
    }
    ofstream fout("output.txt");
    int testCases;
    if (!(fin >> testCases)) return 0;
    DCEL dcel;

    while(testCases--) {
        dcel.clear();
        g_dcel=&dcel;
        g_points=&dcel.points;

        int n;
        if (!(fin >> n)) break;

        vector<pair<double,double>> tmp_outer(n);
        for (int i=0;i<n;i++) fin >> tmp_outer[i].first >> tmp_outer[i].second;

        if (getSignedArea(tmp_outer)<0) reverse(tmp_outer.begin(), tmp_outer.end());

        vector<point*> outer(n);
        for (int i=0;i<n;i++) outer[i]=dcel.join_pt(tmp_outer[i].first, tmp_outer[i].second);
        dcel.add_ring(outer);

        int numHoles;fin >> numHoles;
        dcel.h_count=numHoles;
        for (int hi=0;hi<numHoles;hi++) {
            int k;fin >> k;
            vector<pair<double,double>> tmp_hole(k);
            for (int i=0;i<k;i++) fin >>tmp_hole[i].first >> tmp_hole[i].second;

            if (getSignedArea(tmp_hole)>0) reverse(tmp_hole.begin(), tmp_hole.end());

            vector<point*> hole(k);
            for (int i=0;i<k;i++) {
                hole[i]=dcel.join_pt(tmp_hole[i].first, tmp_hole[i].second);
                hole[i]->is_hole=true;
            }
            dcel.add_ring(hole);
        }

        int n_original_total=(int)dcel.points.size();
        dcel.n_original=n_original_total;

        dcel.orig_prev.assign(n_original_total, -1);
        dcel.orig_next.assign(n_original_total, -1);
        for(auto& edge : dcel.original_edges) {
            dcel.orig_next[edge.first]=edge.second;
            dcel.orig_prev[edge.second]=edge.first;
        }

        make_monotone(n_original_total);
        auto mono_polys=get_monotone_polygon
    (dcel);
        for (auto& poly : mono_polys) triangulate_monotone(dcel, poly);

        vector<int> color=dcel.three_color();
        vector<int> guards[3];
        
        for (int i=0;i<n_original_total;i++) {
            if (color[i] >=0 && color[i]<3) guards[color[i]].push_back(i);
        }
        
        int min_color=0;
        if (guards[1].size()<guards[min_color].size()) min_color=1;
        if (guards[2].size()<guards[min_color].size()) min_color=2;

        for (auto& [a, b, c] : dcel.final_triangles) {
            fout << dcel.points[a]->x << " " << dcel.points[a]->y << " "
                 << dcel.points[b]->x << " " << dcel.points[b]->y << " "
                 << dcel.points[c]->x << " " << dcel.points[c]->y << "\n";
        }
        
        fout << "GUARDS\n";
        fout << guards[min_color].size() << "\n";
        for(int idx : guards[min_color]) {
             fout << dcel.points[idx]->x << " " << dcel.points[idx]->y << "\n";
        }
        fout << "END_TC\n";
    }
    return 0;
}