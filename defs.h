#include<bits/stdc++.h>
using namespace std;

struct point
{
    double x,y;
    point* prev;
    point* nxt;
};


vector<point> outer_poly;
vector<vector<point>> inner_hole;