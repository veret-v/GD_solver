#ifndef POINT
#define POINT

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>

struct Point
{
    double x;
    double y;

    Point(const double _x, const double _y) : x(_x), y(_y) {};
    Point() : x(0.), y(0.) {};

    Point operator+(const Point _other) const { return Point(x + _other.x, y + _other.y); }
    Point operator/(const double _val) const { return Point(x / 2, y / 2); }
    Point operator+=(const Point _other) const { return operator+(_other); }

    void show();
};
#endif