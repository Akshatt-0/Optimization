#ifndef GEOMETRY_UTILS_H
#define GEOMETRY_UTILS_H

#include <cmath>

struct Point {
    double x, y;
    Point(double _x, double _y) : x(_x), y(_y) {}
};

double distance(const Point& a, const Point& b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

#endif // GEOMETRY_UTILS_H