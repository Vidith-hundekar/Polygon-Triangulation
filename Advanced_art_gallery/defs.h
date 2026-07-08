#ifndef DEFS_H
#define DEFS_H

#include<bits/stdc++.h>
struct half_ed;
struct face;

struct point {
    double x, y;
    half_ed* rand_one_ed;
    bool is_hole;           
    bool is_rightmost;      
    int idx;                

    point(double x_i = 0, double y_i = 0) {
        x = x_i; y = y_i;
        rand_one_ed = nullptr;
        is_hole = false; is_rightmost = false;
        idx = -1;
    }
};

struct half_ed {
    point* st;
    half_ed* next;
    half_ed* prev;
    half_ed* reverse_ed;
    face* incident_face;
    half_ed() {
        st = nullptr; next = nullptr; prev = nullptr;
        reverse_ed = nullptr; incident_face = nullptr;
    }
};

struct face {
    half_ed* rand_one_ed;
    face() { rand_one_ed = nullptr; }
};

#endif
