#pragma once

#include "Point.h"

#include "MathTools.h"

#include <iostream>
#include <vector>


class Cube{
private:
    Point bottomLeft1;
    Point bottomRight1;
    Point topLeft1;
    Point topRight1;
    Point bottomLeft2;
    Point bottomRight2;
    Point topLeft2;
    Point topRight2;

public:
    Cube() {}
    Cube(const Point& bl1, const Point& br1, const Point& tl1, const Point& tr1, const Point& bl2, const Point& br2, const Point& tl2, const Point& tr2)
            : bottomLeft1(bl1), bottomRight1(br1), topLeft1(tl1), topRight1(tr1), bottomLeft2(bl2), bottomRight2(br2), topLeft2(tl2), topRight2(tr2) {}
    void printVertices() const {
        std::cout << "Bottom Left 1: ";
        bottomLeft1.print();
        std::cout << std::endl;

        std::cout << "Bottom Right 1: ";
        bottomRight1.print();
        std::cout << std::endl;

        std::cout << "Top Left 1: ";
        topLeft1.print();
        std::cout << std::endl;

        std::cout << "Top Right 1: ";
        topRight1.print();
        std::cout << std::endl;

        std::cout << "Bottom Left 2: ";
        bottomLeft2.print();
        std::cout << std::endl;

        std::cout << "Bottom Right 2: ";
        bottomRight2.print();
        std::cout << std::endl;

        std::cout << "Top Left 2: ";
        topLeft2.print();
        std::cout << std::endl;

        std::cout << "Top Right 2: ";
        topRight2.print();
        std::cout << std::endl;
    }
    Point findBotomLeft(){
        std::vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};
        Point currentBottomLeft = points[0];

        for (int i = 1; i < points.size(); ++i) {
            bool isXClose = isDifferenceLessThan(points[i].x, currentBottomLeft.x, 7);
            bool isYClose = isDifferenceLessThan(points[i].y, currentBottomLeft.y, 7);
            bool isZClose = isDifferenceLessThan(points[i].z, currentBottomLeft.z, 7);

            if ((isXClose && isZClose && points[i].y < currentBottomLeft.y) ||
                (isYClose && isZClose && points[i].x < currentBottomLeft.x) ||
                (isXClose && isYClose && points[i].z < currentBottomLeft.z) ||
                (isXClose && points[i].z < currentBottomLeft.z && points[i].y < currentBottomLeft.y) ||
                (isYClose && points[i].z < currentBottomLeft.z && points[i].x < currentBottomLeft.x) ||
                (isZClose && points[i].y < currentBottomLeft.y && points[i].x < currentBottomLeft.x) ||
                (points[i].y < currentBottomLeft.y && points[i].x < currentBottomLeft.x && points[i].z < currentBottomLeft.z))
            {
                currentBottomLeft = points[i];
            }}
        return currentBottomLeft;
    }
    Point findCenter(){
        std::vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};
        double xcoordinate=0;
        double ycoordinate=0;
        double zcoordinate=0;
        for(int i=0;i<8;i++){
            xcoordinate +=points[i].x;
            ycoordinate +=points[i].y;
            zcoordinate +=points[i].z;
        }
        return Point{xcoordinate/8, ycoordinate/8, zcoordinate/8};
    }
    Point findHeight(){
        std::vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};
        Point currentTop = points[0];

        for (int i = 1; i < points.size(); ++i) {

            if (( points[i].y < currentTop.y))
            {
                currentTop = points[i];
            }}
        return currentTop;
    }
    void findMinMax(Point& min, Point& max) {
        std::vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};

        min.x = min.y = min.z = std::numeric_limits<double>::max();
        max.x = max.y = max.z = std::numeric_limits<double>::lowest();

        for (const auto& vertex : points) {
            if (vertex.x < min.x) min.x = vertex.x;
            if (vertex.y < min.y) min.y = vertex.y;
            if (vertex.z < min.z) min.z = vertex.z;

            if (vertex.x > max.x) max.x = vertex.x;
            if (vertex.y > max.y) max.y = vertex.y;
            if (vertex.z > max.z) max.z = vertex.z;
        }
    }
    std::vector<double> getAngleStretch(Cube cub) {
        Point min, max;
        cub.findMinMax(min, max);
        std::vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};

        double width = max.x - min.x;
        double height = max.y - min.y;
        double depth = max.z - min.z;

        double minDimension = std::min({width, height, depth});

        return {minDimension/width, minDimension/height, minDimension/depth};
    }
    void translate(const Point v) {
        bottomLeft1 += static_cast<Point>(v);
        bottomRight1 += static_cast<Point>(v);
        topLeft1 += static_cast<Point>(v);
        topRight1 += static_cast<Point>(v);
        bottomLeft2 += static_cast<Point>(v);
        bottomRight2 += static_cast<Point>(v);
        topLeft2 += static_cast<Point>(v);
        topRight2 += static_cast<Point>(v);
    }
    void angle_translate(const Point v) {
        bottomLeft1 *= static_cast<Point>(v);
        bottomRight1 *= static_cast<Point>(v);
        topLeft1 *= static_cast<Point>(v);
        topRight1 *= static_cast<Point>(v);
        bottomLeft2 *= static_cast<Point>(v);
        bottomRight2 *= static_cast<Point>(v);
        topLeft2 *= static_cast<Point>(v);
        topRight2 *= static_cast<Point>(v);

    }
    void checkVerticesAndSetZero(int n) {
        bottomLeft1.checkAndSetZero(n);
        bottomRight1.checkAndSetZero(n);
        topLeft1.checkAndSetZero(n);
        topRight1.checkAndSetZero(n);
        bottomLeft2.checkAndSetZero(n);
        bottomRight2.checkAndSetZero(n);
        topLeft2.checkAndSetZero(n);
        topRight2.checkAndSetZero(n);
    }


};