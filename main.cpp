#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <ctime>
#include <fstream>
#include "SFML/Graphics/CircleShape.hpp"
#include <SFML/Graphics.hpp>
#include <sstream>
#include <iomanip>

using namespace std;

struct Point {
    double x, y, z;
    void print() const {
        cout << "(" << x << "," << y << "," << z << ")";
    }
    Point& operator+=(const Point& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    Point& operator*=(const Point& other) {
        x *= other.x;
        y *= other.y;
        z *= other.z;
        return *this;
    }
    void checkAndSetZero(int n);


};
class PlaneEquation {
public:
    double a, b, c, d;

    PlaneEquation(double a, double b, double c, double d) : a(a), b(b), c(c), d(d) {}

    static PlaneEquation buildPlaneEquation(double x1, double y1, double z1, double x2, double y2, double z2) {
        double a = y1 - y2;
        double b = x2 - x1;
        double c = 0;
        double d = -(a * x1 + b * y1);
        return PlaneEquation(a, b, c, d);
    }
    vector<double> getPlaneParams(){return {a,b,c,d};}
};
class Sphere {
private:
    double r;
    vector<Point> p;
    Point c;

public:
    Sphere(double radius, vector<Point> points, Point center)
            : r(radius), p(points), c(center) {}
    double getRadius() const { return r; }
    vector<Point> getPoints() const { return p; }
    Point getCenter() const { return c; }
};
class Circle {
private:
    double r;
    vector<Point> p;
    Point c;

public:
    Circle(double radius, vector<Point> points, Point center)
            : r(radius), p(points), c(center) {}
    double getRadius() const { return r; }
    vector<Point> getPoints() const { return p; }
    Point getCenter() const { return c; }
};
bool isDifferenceLessThan(double a, double b, int n) {
    double threshold = std::pow(10, -n);
    return std::abs(a - b) < threshold;
}

class LineEquation {
public:
    double a, b, c, a0, b0, c0;

    LineEquation(double a, double b, double c, double a0, double b0, double c0) : a(a), b(b), c(c), a0(a0), b0(b0), c0(c0) {}
    //x-a0/a=y-b0/b
    //x2 = a-a0
    string getCanonicalEquation() const {
        string eq = "(x - " + to_string(a0) + ") / " + to_string(a) +
                    " = (y - " + to_string(b0) + ") / " + to_string(b) +
                    " = (z - " + to_string(c0) + ") / " + to_string(c);
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
    Cube(Point& bl1, Point& br1, Point& tl1, Point& tr1, Point& bl2, Point& br2, Point& tl2, Point& tr2)
    : bottomLeft1(bl1), bottomRight1(br1), topLeft1(tl1), topRight1(tr1), bottomLeft2(bl2), bottomRight2(br2), topLeft2(tl2), topRight2(tr2) {}
    void printVertices() const {
        cout << "Bottom Left 1: ";
        bottomLeft1.print();
        cout << endl;

        cout << "Bottom Right 1: ";
        bottomRight1.print();
        cout << endl;

        cout << "Top Left 1: ";
        topLeft1.print();
        cout << endl;

        cout << "Top Right 1: ";
        topRight1.print();
        cout << endl;

        cout << "Bottom Left 2: ";
        bottomLeft2.print();
        cout << endl;

        cout << "Bottom Right 2: ";
        bottomRight2.print();
        cout << endl;

        cout << "Top Left 2: ";
        topLeft2.print();
        cout << endl;

        cout << "Top Right 2: ";
        topRight2.print();
        cout << endl;
    }
    Point findBotomLeft(){
        vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};
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
        vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};
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
        vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};
        Point currentTop = points[0];

        for (int i = 1; i < points.size(); ++i) {

            if (( points[i].y < currentTop.y))
            {
                currentTop = points[i];
            }}
        return currentTop;
    }
    void findMinMax(Point& min, Point& max) {
        vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};

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
    vector<double> getAngleStretch(Cube cub) {
        Point min, max;
        cub.findMinMax(min, max);
        vector<Point> points= {bottomLeft1, bottomRight1, topLeft1, topRight1, bottomLeft2, bottomRight2, topLeft2, topRight2};

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
void Point::checkAndSetZero(int n) {
    if (isDifferenceLessThan(x, 0.0, n)) x = 0.0;
    if (isDifferenceLessThan(y, 0.0, n)) y = 0.0;
    if (isDifferenceLessThan(z, 0.0, n)) z = 0.0;
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



int main() {
    sf::RenderWindow window(sf::VideoMode(1000, 1000), "15");
    window.setFramerateLimit(60);
    sf::View view;
    view.setCenter(0, 0);
    view.setSize(40, 40);

    vector<vector<double>> all_frequency_marker;

    int n = 6;
    std::vector<int> ms{100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    vector<Point> points;
    for (int i = 0; i < n; i++) {
        double x = randNorm(0.0, 1.0);
        double y = randNorm(0.0, 1.0);
        double z = randNorm(0.0, 1.0);
        points.push_back({x, y, z});

    }
    //points = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {3, 6, 2}, {5, 1, 8}, {9, 4, 7}};
    //1.
    Point maxDistPoint1;
    Point maxDistPoint2;
    double maxDist = 0;

    pair<int, int> max_length_point_place;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = points[i].x - points[j].x;
            double dy = points[i].y - points[j].y;
            double dz = points[i].z - points[j].z;
            double dist = sqrt(dx * dx + dy * dy + dz * dz);

            if (dist > maxDist) {
                maxDist = dist;
                maxDistPoint1 = points[i];
                maxDistPoint2 = points[j];
                max_length_point_place = {i, j};
            }
        }
    }

    //2.

    //3
    LineEquation L = LineEquation::buildLineEquation(points[max_length_point_place.first], points[max_length_point_place.second]);

    //4.
    /*LineEquation L1 = L.parallelLine(xr);
    LineEquation L2 = L.parallelLine(xq);*/
    //6.

    cout<<"original?"<<endl;
    points[max_length_point_place.first].print();cout<<endl;
    points[max_length_point_place.second].print();cout<<endl;
    cout<<endl;

    //max_width_point.first.print();cout << endl;
    //max_width_point.second.print();cout << endl;
    cout<<L.getCanonicalEquation()<<endl;
    /*cout<<L1.getCanonicalEquation()<<endl;
    cout<<L2.getCanonicalEquation()<<endl;*/

    double pi = atan(1) * 4;
    vector<double> angles;
    for(int i=0;i<3;i++) angles.push_back( L.find_angle(i));
    cout<<"angles: "<<angles[0]<<", "<<angles[1]<<", "<<angles[2]<<endl;
    //7.
    //навколо z



    for (int i = 0; i < n; i++) {
        double x = points[i].x;
        double y = points[i].y;
        double z = points[i].z;
        points[i].x = (x * cos(-angles[2])) - (y * sin(-angles[2]));
        points[i].y = (x * sin(-angles[2])) + y * cos(-angles[2]);
        points[i].z = z;
    }

    L = LineEquation::buildLineEquation(points[max_length_point_place.first], points[max_length_point_place.second]);
    cout<<L.getCanonicalEquation()<<endl;

    for(int i=0;i<3;i++) angles[i] = ( L.find_angle(i));
    cout<<"angles: "<<angles[0]<<", "<<angles[1]<<", "<<angles[2]<<endl;

    for (int i = 0; i < n; i++) {
        double x = points[i].x;
        double y = points[i].y;
        double z = points[i].z;
        points[i].x = (x * cos(  -angles[1])) + (z * sin(  -angles[1]));
        points[i].y = y;
        points[i].z = (-x * sin(  -angles[1])) + (z * cos(  -angles[1]));
    }


    L = LineEquation::buildLineEquation(points[max_length_point_place.first], points[max_length_point_place.second]);
    auto furthestPair = findFurthestPair(points, L);
    Point xr = furthestPair.first;
    Point xq = furthestPair.second;
    pair<Point, Point> max_width_point = {xr, xq};
    cout<<L.getCanonicalEquation()<<endl;

    for(int i=0;i<3;i++) angles[i] = (L.find_angle(i));
    cout<<"angles: "<<angles[0]<<", "<<angles[1]<<", "<<angles[2]<<endl;

    LineEquation L1 = L.parallelLine(max_width_point.first);
    LineEquation L2 = L.parallelLine(max_width_point.second);
    //L1.findIntersectionWithPerpendicularLine();
    //L3 = L3.find_rotated_line(L3, angle);
    //L4 = L4.find_rotated_line(L4, angle);
    cout<<"original?"<<endl;
    points[max_length_point_place.first].print();cout<<endl;
    points[max_length_point_place.second].print();cout<<endl;
    cout<<endl;
    max_width_point.first.print();cout << endl;
    max_width_point.second.print();cout << endl;
   cout<<L1.getCanonicalEquation()<<endl;
   cout<<L2.getCanonicalEquation()<<endl;
    //LineEquation L3 = L.perpendicularLine(maxDistPoint1);
    //LineEquation L4 = L.perpendicularLine(maxDistPoint2);

    //6.
    vector<Point> rectangular_points = {L1.findIntersectionWithPerpendicularLine(points[max_length_point_place.first]),
                                        L2.findIntersectionWithPerpendicularLine(points[max_length_point_place.first]),
                                        L1.findIntersectionWithPerpendicularLine(points[max_length_point_place.second]),
                                        L2.findIntersectionWithPerpendicularLine(points[max_length_point_place.second])};
    Cube cub{rectangular_points[0],
             rectangular_points[1],
             rectangular_points[2],
             rectangular_points[3],
             get_all_corners(points[max_length_point_place.first], make_pair(rectangular_points[0], rectangular_points[1])).first,
             get_all_corners(points[max_length_point_place.first], make_pair(rectangular_points[0], rectangular_points[1])).second,

             get_all_corners(points[max_length_point_place.second], make_pair(rectangular_points[2], rectangular_points[3])).first,
             get_all_corners(points[max_length_point_place.second], make_pair(rectangular_points[2], rectangular_points[3])).second};


    cub.printVertices();

    Rectangle rect{};
    Point transfer = {0 - cub.findBotomLeft().x, 0 - cub.findBotomLeft().y, 0 - cub.findBotomLeft().z};

    for (int i = 0; i < n; i++) {
        points[i] += transfer;
    }
    max_width_point.first += transfer;
    max_width_point.second += transfer;

    cub.translate(transfer);
    cub.printVertices();
    vector<double> alpha = cub.getAngleStretch(cub);
    for (int i = 0; i < n; i++) {
        points[i] *= {alpha[0], alpha[1], alpha[2]};
    }
    max_width_point.first *= {alpha[0], alpha[1], alpha[2]};
    max_width_point.second *= {alpha[0], alpha[1], alpha[2]};


    cub.angle_translate({alpha[0], alpha[1], alpha[2]});
    cub.checkVerticesAndSetZero(7);
    cub.printVertices();

    Point center = cub.findCenter();

    vector<double> distances(n);
    for (int i = 0; i < points.size(); i++) {
        distances[i] = sqrt(pow(center.x - points[i].x, 2) + pow(center.y - points[i].y, 2)+ pow(center.z - points[i].z, 2));
    }

    double radius = *std::max_element(distances.begin(), distances.end());
    size_t max_index = std::distance(distances.begin(), std::max_element(distances.begin(), distances.end()));

    vector<Point> point_previous = points;
    points[3].print();cout<<endl;
    Sphere sphere(radius, points, center);

    sf::CircleShape circleDrawing(radius);
    sf::CircleShape circleDrawingUnshuffled(radius);
    circleDrawing.setScale((1 / alpha[0]), 1);

    for (int i = 0; i < points.size(); i++) {
        sphere.getPoints()[i].x *= 1 / alpha[0];
        points[i].x *= 1 / alpha[0];
        if (points[i].x < 0 || points[i].y < 0) {
            cout << points[i].x << " " << points[i].y;
        }
    }

    double a = sqrt(pow(sphere.getPoints()[max_width_point.first.x].x - center.x, 2) +
                    pow(sphere.getPoints()[max_width_point.first.x].y - center.y, 2) +
                    pow(sphere.getPoints()[max_width_point.first.x].z - center.z, 2));
    double b = sqrt(pow(sphere.getPoints()[max_length_point_place.first].x - center.x, 2) +
                    pow(sphere.getPoints()[max_length_point_place.first].y - center.y, 2)+
                            pow(sphere.getPoints()[max_width_point.first.x].z - center.z, 2));
    double c = sqrt(pow(sphere.getPoints()[max_length_point_place.first].x - center.x, 2) +
                    pow(sphere.getPoints()[max_length_point_place.first].y - center.y, 2)+
                            pow(sphere.getPoints()[max_length_point_place.first].z - center.z, 2));
    cout << "ellipse equation = (x - " << center.x << ")**2 / " << pow(b, 2) << " + " << " (y - "
         << center.y << ")**2 / " << pow(a, 2) << " = 1" << endl;

    vector<double> ms_frequency_marker;
    for (auto &&m : ms) {
        std::vector<Point> m_points;
        for (int i = 0; i < m; i++) {
            double x = randNorm(0.0, 1.0);
            double y = randNorm(0.0, 1.0);
            m_points.push_back({x, y});
        }

        int number_points_in_ellipse = 0;
        double frequency = frequencyCheck(sphere, m_points, a, b);
        ms_frequency_marker.push_back(frequency);
        cout << frequency << " ";
    }
    all_frequency_marker.push_back(ms_frequency_marker);
    cout << endl;

    ofstream outfile("results.txt");
    if (outfile.is_open()) {
        for (int i = 0; i < all_frequency_marker.size(); i++) {
            for (int j : all_frequency_marker[i]) {
                outfile << j << " ";
            }
            outfile << std::endl;
        }
        outfile.close();
    }

    while (window.isOpen()) {
        sf::Event event;
        window.setView(view);

        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        window.clear(sf::Color::Black);

        sf::Vertex centerDrawing(sf::Vector2f(center.x, center.y), sf::Color::Red);

        circleDrawing.setPosition((center.x - radius / alpha[0]), (center.y - radius));
        circleDrawing.setFillColor(sf::Color::Transparent);
        circleDrawing.setOutlineThickness(0.1);
        circleDrawing.setOutlineColor(sf::Color(250, 150, 100));
        circleDrawingUnshuffled.setPosition(center.x - radius, center.y - radius);
        circleDrawingUnshuffled.setOutlineColor(sf::Color::Blue);
        circleDrawingUnshuffled.setFillColor(sf::Color::Transparent);
        circleDrawingUnshuffled.setOutlineThickness(0.1);

        sf::VertexArray axes(sf::Lines);

        axes.append(sf::Vertex(sf::Vector2f(-800, 0), sf::Color::Red));
        axes.append(sf::Vertex(sf::Vector2f(800, 0), sf::Color::Red));
        axes.append(sf::Vertex(sf::Vector2f(0, -800), sf::Color::Red));
        axes.append(sf::Vertex(sf::Vector2f(0, 800), sf::Color::Red));

        axes.append(sf::Vertex(sf::Vector2f(1, -1), sf::Color::Red));
        axes.append(sf::Vertex(sf::Vector2f(1, 1), sf::Color::Red));
        for (size_t i = 0; i < points.size(); ++i) {
            sf::CircleShape pointShape(0.05);
            sf::CircleShape point2shape(0.05);
            point2shape.setFillColor(sf::Color::Red);
            pointShape.setFillColor(sf::Color::White);
            pointShape.setPosition(point_previous[i].x, point_previous[i].y);
            point2shape.setPosition(points[i].x, points[i].y);
            window.draw(point2shape);
            window.draw(pointShape);
        }

        window.draw(axes);
        window.draw(&centerDrawing, 50, sf::Points);
        window.draw(circleDrawing);
        window.draw(circleDrawingUnshuffled);

        window.display();
    }

    return 0;
}

