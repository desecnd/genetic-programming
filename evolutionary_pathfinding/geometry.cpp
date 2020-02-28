/*
* source: Geometry for Competitve Programming guide bu vlecomte on codeforces
*/
#include "geometry.h"

using geo::Point;
using geo::Line;
using geo::Circle;
using geo::T;

inline double geo::PI() { return std::atan(1)*4; }

bool geo::operator==(Point a, Point b) { return a.x == b.x && a.y == b.y; }
bool geo::operator!=(Point a, Point b) { return !(a == b); }

std::ostream& geo::operator<<(std::ostream& out, Point p) {
    return out << "(" << p.x << "," << p.y << ")";
}

Point geo::perp(geo::Point p) { return { -p.y, p.x }; }
T geo::sq(Point p) { return p.x*p.x + p.y*p.y; }
double geo::abs(Point p) { return std::sqrt(geo::sq(p)); } 

T geo::dot(Point v, Point w) { return v.x * w.x + v.y * w.y; }
double geo::angle(Point v, Point w) {
    double cosTheta = geo::dot(v, w) / geo::abs(v) / geo::abs(w);
    return std::acos( std::max(-1.0, std::min(1.0, cosTheta))); 
}

T geo::cross(Point v, Point w) { return v.x * w.y - w.x * v.y; }
T geo::orient(Point a, Point b, Point c) { return geo::cross(b - a, c - a); }
double geo::orientedAngle(Point a, Point b, Point c) {
    if ( geo::orient(a, b, c) >= 0 ) return geo::angle(b - a, c - a);
    else return 2 * geo::PI() - geo::angle(b - a, c - a); 
}
bool geo::inAngle(Point a, Point b, Point c, Point p) {
    // std::assert( geo::orient(a,b,c) != 0 ); // no angle
    if ( geo::orient(a, b, c) < 0 ) std::swap(b, c);
    return geo::orient(a, b, p) >= 0 && geo::orient(a, c, p) <= 0;
}


bool geo::inter(Line a, Line b, Point& out) {
    if ( geo::cross(a.v, b.v) == 0 ) return false;
    out = (b.v * a.c - a.v * b.c) / geo::cross(a.v, b.v);
    return true;
} 

double geo::segPoint(Point a, Point b, Point p) {
    if ( a != b ) {
        Line l{a, b};
        if ( l.cmpProj(a, p) && l.cmpProj(p, b) ) return l.dist(p);
    }
    return std::min(abs(p - a), abs(p - b));
}

