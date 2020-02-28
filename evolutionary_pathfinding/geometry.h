/*
* source: Geometry for Competitve Programming guide bu vlecomte on codeforces
*/
#ifndef GEOMETRY_8656167_H
#define GEOMETRY_8656167_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>

namespace geo { 
using T = long double;
inline double PI();

struct Point {
    T x, y;
    Point(T x = 0, T y = 0) : x(x), y(y) {}

    Point operator+(Point p) { return { x + p.x, y + p.y }; }
    Point operator-(Point p) { return { x - p.x, y - p.y }; }
    Point operator*(T val) { return { x * val, y * val }; }
    Point operator/(T val) { return { x / val, y / val }; }
};

bool operator==(Point a, Point b); 
bool operator!=(Point a, Point b);

std::ostream& operator<<(std::ostream& out, Point p);

Point perp(Point p);
T sq(Point p);
double abs(Point p);

T dot(Point v, Point w);
double angle(Point v, Point w);

T cross(Point v, Point w);
T orient(Point a, Point b, Point c);
double orientedAngle(Point a, Point b, Point c);
bool inAngle(Point a, Point b, Point c, Point p);

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

bool inter(Line a, Line b, Point& out);
double segPoint(Point a, Point b, Point p);

struct Circle {
    Point o; double r;
    Circle() : o(), r(0) {}
    Circle(Point p, double r) : o(p), r(r) {}
    Circle(Point a, Point b, Point c) {
        Line perpAB{ Line(a, b).perpThrough((a + b) / 2) };
        Line perpAC{ Line(a, c).perpThrough((a + c) / 2) };
        // std::assert(inter(perpAB, perpAC, o));
        r = abs(a - o);
    }
};
};

#endif