#pragma once

#include "Point.h"

#include <vector>

class Sphere {
private:
    double r;
    std::vector<Point> p;
    Point c;

public:
    Sphere(double radius, std::vector<Point> points, Point center)
            : r(radius), p(points), c(center) {}
    double getRadius() const { return r; }
    std::vector<Point> getPoints() const { return p; }
    Point getCenter() const { return c; }
};
