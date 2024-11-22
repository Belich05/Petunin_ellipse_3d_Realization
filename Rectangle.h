#pragma once
#include <iostream>
#include "Point.h"
#include <vector>
using namespace::std;

class Rectangle {
private:
    Point bottomLeft;
    Point bottomRight;
    Point topLeft;
    Point topRight;

public:
    Rectangle() {}
    Rectangle(Point& bl, Point& br, Point& tl, Point& tr)
            : bottomLeft(bl), bottomRight(br), topLeft(tl), topRight(tr) {}

    Point getBottomLeft() const { return bottomLeft; }
    Point getBottomRight() const { return bottomRight; }
    Point getTopLeft() const { return topLeft; }
    Point getTopRight() const { return topRight; }
    void setPoints(double alpha, Point transfer) {
        Point the_most_bottomLeft = findBottomLeft();
        Point move_vector = transfer;
        bottomLeft.x += move_vector.x; bottomRight.x +=move_vector.x; topLeft.x += move_vector.x; topRight.x += move_vector.x;
        bottomLeft.y += move_vector.y; bottomRight.y +=move_vector.y; topLeft.y += move_vector.y; topRight.y += move_vector.y;
    }
    Point findBottomLeft() {
        std::vector<Point> points = {bottomLeft, bottomRight, topLeft, topRight};
        Point currentBottomLeft =points[0];
        for (int i=0; i<points.size(); i++) {
            if((isDifferenceLessThan(points[i].x, currentBottomLeft.x, 7) && (points[i].y < currentBottomLeft.y))||(points[i].x < currentBottomLeft.x) ) {
                currentBottomLeft = points[i];
            }
        }
        return currentBottomLeft;
    }
    Point findTopRight() {
        std::vector<Point> points = {bottomLeft, bottomRight, topLeft, topRight};
        Point currentBottomLeft =points[0];
        for (int i=0; i<points.size(); i++) {
            if(!(isDifferenceLessThan(points[i].x, 0, 7)) && !(isDifferenceLessThan(points[i].y, 0, 7))) currentBottomLeft = points[i];
        }
        return currentBottomLeft;
    }

    double getWidth() const {
        return sqrt(pow(bottomRight.x - bottomLeft.x, 2) +
                    pow(bottomRight.y - bottomLeft.y, 2));
    }

    double getHeight() const {
        return sqrt(pow(topLeft.x - bottomLeft.x, 2) +
                    pow(topLeft.y - bottomLeft.y, 2));
    }
};