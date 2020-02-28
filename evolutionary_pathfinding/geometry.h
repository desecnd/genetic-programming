/*
* source: Geometry for Competitve Programming guide bu vlecomte on codeforces
*/

#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <tuple>
using namespace std;

using T = long double;
inline double PI() { return atan(1)*4; }
struct Point {
    T x, y;
    Point(T x = 0, T y = 0) : x(x), y(y) {}

    Point operator+(Point p) { return { x + p.x, y + p.y }; }
    Point operator-(Point p) { return { x - p.x, y - p.y }; }
    Point operator*(T val) { return { x * val, y * val }; }
    Point operator/(T val) { return { x / val, y / val }; }
};

bool operator==(Point a, Point b) { return a.x == b.x && a.y == b.y; }
bool operator!=(Point a, Point b) { return !(a == b); }

ostream& operator<<(ostream& out, Point p) {
    return out << "(" << p.x << "," << p.y << ")";
}

Point perp(Point p) { return { -p.y, p.x }; }
T sq(Point p) { return p.x*p.x + p.y*p.y; }
double abs(Point p) { return sqrt(sq(p)); } 

T dot(Point v, Point w) { return v.x * w.x + v.y * w.y; }
double angle(Point v, Point w) {
    double cosTheta = dot(v, w) / abs(v) / abs(w);
    return acos( max(-1.0, min(1.0, cosTheta))); 
}

T cross(Point v, Point w) { return v.x * w.y - w.x * v.y; }
T orient(Point a, Point b, Point c) { return cross(b - a, c - a); }
double orientedAngle(Point a, Point b, Point c) {
    if ( orient(a, b, c) >= 0 ) return angle(b - a, c - a);
    else return 2 * PI() - angle(b - a, c - a); 
}
bool inAngle(Point a, Point b, Point c, Point p) {
    assert( orient(a,b,c) != 0 ); // no angle
    if ( orient(a, b, c) < 0 ) swap(b, c);
    return orient(a, b, p) >= 0 && orient(a, c, p) <= 0;
}

struct Line {
    Point v; T c;
    Line(Point v, T c) : v(v), c(c) {} 
    Line(T a, T b, T c) : v(b, -a), c(c) {} 
    Line(Point p, Point q) : v(q - p), c(cross(v, p)) {} 

    T side(Point p) { return cross(v, p) - c; }
    Line perpThrough(Point p) { return { p , p + perp(v) }; }
    bool cmpProj(Point a, Point b) { return dot(v, a) < dot(v, b); }
    double dist(Point p) { return abs(side(p)) / abs(v); }
};

bool inter(Line a, Line b, Point& out) {
    if ( cross(a.v, b.v) == 0 ) return false;
    out = (b.v * a.c - a.v * b.c) / cross(a.v, b.v);
    return true;
} 

double segPoint(Point a, Point b, Point p) {
    if ( a != b ) {
        Line l{a, b};
        if ( l.cmpProj(a, p) && l.cmpProj(p, b) ) return l.dist(p);
    }
    return min(abs(p - a), abs(p - b));
}

struct Circle {
    Point o; double r;
    Circle() : o(), r(0) {}
    Circle(Point p, double r) : o(p), r(r) {}
    Circle(Point a, Point b, Point c) {
        Line perpAB{ Line(a, b).perpThrough((a + b) / 2) };
        Line perpAC{ Line(a, c).perpThrough((a + c) / 2) };
        assert(inter(perpAB, perpAC, o));
        r = abs(a - o);
    }
};