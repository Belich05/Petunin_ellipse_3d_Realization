#pragma once
#include<string>
#include<iostream>

#include "Point.h"

class LineEquation {
public:
    double a, b, c, a0, b0, c0;

    LineEquation(double a, double b, double c, double a0, double b0, double c0) : a(a), b(b), c(c), a0(a0), b0(b0), c0(c0) {}
    //x-a0/a=y-b0/b
    //x2 = a-a0
    std::string getCanonicalEquation() const {
        std::string eq = "(x - " + std::to_string(a0) + ") / " + std::to_string(a) +
                    " = (y - " + std::to_string(b0) + ") / " + std::to_string(b) +
                    " = (z - " + std::to_string(c0) + ") / " + std::to_string(c);
        return eq;
    }


    static LineEquation buildLineEquation(Point a1, Point a2) {
        double x1 = a1.x;
        double y1 = a1.y;
        double z1 = a1.z;
        double x2 = a2.x;
        double y2=  a2.y;
        double z2 = a2.z;

        double a0 = x1;
        double b0 = y1;
        double c0 = z1;
        double a = x2-  x1;
        double b = y2-  y1;
        double c = z2-  z1;
        return LineEquation(a,b,c, a0,b0,c0);
    }

    Point crossProduct(Point v1, Point v2) {
        Point result;
        result.x = v1.y * v2.z - v1.z * v2.y;
        result.y = v1.z * v2.x - v1.x * v2.z;
        result.z = v1.x * v2.y - v1.y * v2.x;
        return result;
    }

    PlaneEquation planeEquation(Point point) {
        Point vec;
        vec.x = 0;
        vec.y = point.y - b0;
        vec.z = 0 - c0;

        Point normal = crossProduct({a,b,c}, vec);

        // Рівняння площини: A(x - x0) + B(y - y0) + C(z - z0) = 0
        double A = normal.x;
        double B = normal.y;
        double C = normal.z;
        double D = -(A * point.x + B * point.y + C * point.z);


        return PlaneEquation{A,B,C,D};

    }

    double distanceToPointXoZ(Point distant_point) {
        double x0 = distant_point.x;
        double z0 = distant_point.z;

        // Vector v = (x0 - a0, z0 - c0)
        double vx = x0 - a0;
        double vz = z0 - c0;

        // Vector d = (x2 - x1, z2 - z1)
        double dx = a;
        double dz = c;

        // Cross product v x d in XoZ plane (only y-component is considered)
        double cy = vz * dx - vx * dz;

        double crossProductMagnitude = (cy);

        double directionMagnitude = std::sqrt(dx * dx + dz * dz);
        double distance = crossProductMagnitude / directionMagnitude;
        return distance;
    }
    double distanceToPointYoZ(Point distant_point) {
        double y0 = distant_point.y;

        // Vector v = (y0 - b0, z0 - c0)
        double vy = y0 - b0;


        return vy;
    }





    Point findIntersectionWithPerpendicularLine(Point point_through) {
        // Напрямний вектор прямої
        double x1p = point_through.x;
        double y1p = point_through.y;
        double z1p = point_through.z;

        double dx = a;
        double dy = b;
        double dz = c;

        // Точка на прямій
        double x0 = a0;
        double y0 = b0;
        double z0 = c0;

        // Вектор від точки на прямій до точки через яку проходить перпендикуляр
        double vx = x1p - x0;
        double vy = y1p - y0;
        double vz = z1p - z0;

        // Проекція вектора (vx, vy, vz) на напрямний вектор (dx, dy, dz)
        double t = (vx * dx + vy * dy + vz * dz) / (dx * dx + dy * dy + dz * dz);

        // Точка перетину
        double x = x0 + t * dx;
        double y = y0 + t * dy;
        double z = z0 + t * dz;

        std::cout << "Intersection point with the perpendicular line: (" << x << ", " << y << ", " << z << ")" << std::endl;
        return Point{x, y, z};
    }

    LineEquation parallelLine(const Point& point) const {
        double x0 = point.x;
        double y0 =  point.y;
        double z0 = point.z;


        return LineEquation(a, b, c, x0,y0,z0);
    }

    /*LineEquation moveLine(double dx, double dy, double dz) const {
        double newD = d - (a * dx + b * dy + c * dz);
        return LineEquation(a, b, c, d, e, newD);
    }*/

    /*LineEquation find_rotated_line(LineEquation &L, double angle){
        double aba = cos(angle) * a - sin(angle) * b;
        double baba = sin(angle) * a + cos(angle) * b;
        return LineEquation(aba, baba, c, d, e, f);
    }*/

    /*LineEquation find_compressed_line(LineEquation &L, double coef){
        return LineEquation(a, b, c, d, e, f * coef);
    }*/

    double find_angle(int excluded_axis) const {
        //0 = x, 1=y, 2=z
        double acoef = a;
        double bcoef = b;
        double ccoef = c;
        double result;

        if((excluded_axis) == 0) {
            if(isDifferenceLessThan(b, 0, 7)){
                return atan(c/a);
            }
            result = -(atan(c/sqrt(a*a+c*c)));
        }
        else if((excluded_axis) == 1) {
            result = -atan(c/sqrt(a*a+b*b));
        }
        else if((excluded_axis) == 2) {

            result = atan2(b, a);
        }
        if(isDifferenceLessThan(result, 0, 7)) return 0;

        return result;
    }
};