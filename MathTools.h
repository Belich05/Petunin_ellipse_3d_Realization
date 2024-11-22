#pragma once

#include "Point.h"

#include <random>


bool isDifferenceLessThan(double a, double b, int n) {
    double threshold = std::pow(10, -n);
    return std::abs(a - b) < threshold;
}

Point rotateToXAxis(double a, double b, double c) {
    double theta = std::atan2(b, a);
    double phi = std::atan2(c, std::sqrt(a * a + b * b));

    double cosTheta = std::cos(theta);
    double sinTheta = std::sin(theta);
    double cosPhi = std::cos(phi);
    double sinPhi = std::sin(phi);

    double aPrime = cosPhi * cosTheta * a - cosPhi * sinTheta * b + sinPhi * c;
    double bPrime = sinTheta * a + cosTheta * b;
    double cPrime = -sinPhi * cosTheta * a + sinPhi * sinTheta * b + cosPhi * c;

    return {aPrime, bPrime, cPrime};
}


double distanceBetweenPoints(const Point& p1, const Point& p2) {
    return std::sqrt(std::pow(p2.x - p1.x, 2) +
                     std::pow(p2.y - p1.y, 2) +
                     std::pow(p2.z - p1.z, 2));
}

double randNorm(double mean, double stddev) {
    static random_device rd;
    static mt19937 gen(rd());
    static normal_distribution<double> dist(mean, stddev);
    return dist(gen);
}

pair<Point, Point> findFurthestPair(const vector<Point>& points, LineEquation L) {
    pair<Point, Point> furthestPair;
    double maxPositiveDistanceXoZ = 0;
    double maxNegativeDistanceXoZ = 0;
    double maxPositiveDistanceYoZ = 0;
    double maxNegativeDistanceYoZ = 0;

    for (int i = 0; i < points.size(); ++i) {
        auto C = points[i];
        double distanceXoZ = L.distanceToPointXoZ(C);

        if (distanceXoZ > maxPositiveDistanceXoZ && distanceXoZ > 0) {
            maxPositiveDistanceXoZ = distanceXoZ;
            furthestPair.first = C;
        } else if (distanceXoZ < maxNegativeDistanceXoZ && distanceXoZ < 0) {
            maxNegativeDistanceXoZ = distanceXoZ;
            furthestPair.second = C;
        }
    }
    for (int i = 0; i < points.size(); ++i) {
        auto C = points[i];
        double distanceYoZ = - L.distanceToPointYoZ(C);
        if (distanceYoZ > maxPositiveDistanceYoZ && distanceYoZ > 0) {
            maxPositiveDistanceYoZ = distanceYoZ;
            furthestPair.first.y = C.y ;
        } else if (distanceYoZ < maxNegativeDistanceYoZ && distanceYoZ < 0) {
            maxNegativeDistanceYoZ = distanceYoZ;
            furthestPair.second.y = C.y;

        }
    }



    return furthestPair;
}


double frequencyCheck(Sphere& sphere, vector<Point> m, double a, double b){
    int count = 0;
    for (Point val : m) {
        double d = pow(val.x -sphere.getCenter().x,2) /pow(a,2)
                   + pow(val.y -sphere.getCenter().y,2) /pow(b,2);

        if (d <= 1) {
            count++;
        }
    }
    return count;
}
Point get_mirrored_point(Point point, Point dir_vector, bool which_axis_to_mirror) {
    // 1 для y, 0 для z

    if (which_axis_to_mirror == 0) {
        // Відображення відносно осі z
        point.y -= dir_vector.y;
        point.z += dir_vector.z;
    } else {
        // Відображення відносно осі y
        point.y += dir_vector.y;
        point.z -= dir_vector.z;
    }
    return point;
}
pair<Point,Point> get_all_corners(Point insidePoint,  pair<Point,Point> cornerPoints){
    //ми змінюємо знак 1 координати одного з напрямляючих векторів, а потім дивимося відстань до точки яку ми перевірили, вона має бути найменшою
    Point first_dir_vector = {cornerPoints.first.x - insidePoint.x, cornerPoints.first.y - insidePoint.y, cornerPoints.first.z - insidePoint.z};
    Point second_dir_vector = {cornerPoints.second.x - insidePoint.x, cornerPoints.second.y - insidePoint.y, cornerPoints.second.z - insidePoint.z};

    //Point mirrored_first_point = get_mirrored_point(insidePoint,first_dir_vector, 0);
    //    Point mirrored_second_point = get_mirrored_point(insidePoint,second_dir_vector, 1);
    Point mirrored_first_point = {cornerPoints.first.x, cornerPoints.first.y, cornerPoints.second.z};
    Point mirrored_second_point ={cornerPoints.first.x, cornerPoints.second.y, cornerPoints.first.z};
    /*if(cornerPoints.second.x != mirrored_first_point.x){
        mirrored_first_point = get_mirrored_point(insidePoint,first_dir_vector, 1);
        mirrored_second_point = get_mirrored_point(insidePoint,second_dir_vector, 0);
    }*/

    return make_pair(mirrored_first_point,mirrored_second_point);
}
Point rotateVector(double angle, Point initial, Point rotation_vector) { double c = cos(angle);
    double x = initial.x;
    double y = initial.y;
    double z = initial.z;

    double vx = rotation_vector.x;
    double vy = rotation_vector.y;
    double vz = rotation_vector.z;

    double s = sin(angle);

    double length = sqrt(x*x + y*y + z*z);

    x /= length;

    y /= length;

    z /= length;

    double rotationMatrix[3][3] = {

            { c + x*x*(1-c), x*y*(1-c) - z*s, x*z*(1-c) + y*s },

            { y*x*(1-c) + z*s, c + y*y*(1-c), y*z*(1-c) - x*s },

            { z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c) } };

    double newVx = rotationMatrix[0][0]*vx + rotationMatrix[0][1]*vy + rotationMatrix[0][2]*vz;

    double newVy = rotationMatrix[1][0]*vx + rotationMatrix[1][1]*vy + rotationMatrix[1][2]*vz;

    double newVz = rotationMatrix[2][0]*vx + rotationMatrix[2][1]*vy + rotationMatrix[2][2]*vz;
    return Point{newVx, newVy, newVz};
}

