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

#include "Point.h"
#include "PlaneEquation.h"
#include "Sphere.h"
#include "LineEquation.h"
#include "Rectangle.h"
#include "Cube.h"
#include "Circle.h"


using namespace std;



int main() {
    sf::RenderWindow xy_projection_window(sf::VideoMode(1000, 1000), "15");
    xy_projection_window.setFramerateLimit(60);
    sf::View view;
    view.setCenter(0, 0);
    view.setSize(40, 40);

    vector<vector<double>> all_frequency_marker;

    int n = 6;
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
    LineEquation L = LineEquation::buildLineEquation(points[max_length_point_place.first], points[max_length_point_place.second]);


    cout<<"original?"<<endl;
    points[max_length_point_place.first].print();cout<<endl;
    points[max_length_point_place.second].print();cout<<endl;
    cout<<endl;

    cout<<L.getCanonicalEquation()<<endl;

    double pi = atan(1) * 4;
    vector<double> angles;
    for(int i=0;i<3;i++) angles.push_back( L.find_angle(i));
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
    for(int i=0;i<3;i++) angles[i] = (L.find_angle(i));

    LineEquation L1 = L.parallelLine(max_width_point.first);
    LineEquation L2 = L.parallelLine(max_width_point.second);
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
    circleDrawing.setScale((1 / alpha[0]), 1 / alpha[1]);

    for (int i = 0; i < points.size(); i++) {
        sphere.getPoints()[i].x *= 1 / alpha[0];
        points[i].x *= 1 / alpha[0];
        points[i].y *= 1 / alpha[1];
        points[i].z *= 1 / alpha[2];
    }


    while (xy_projection_window.isOpen()) {
        sf::Event event;
        xy_projection_window.setView(view);

        while (xy_projection_window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                xy_projection_window.close();
            }
        }

        xy_projection_window.clear(sf::Color::Black);

        sf::Vertex centerDrawing(sf::Vector2f(center.x, center.y), sf::Color::Red);

        circleDrawing.setPosition((center.x - radius ), (center.y - radius));
        circleDrawing.setFillColor(sf::Color::Transparent);
        circleDrawing.setOutlineThickness(0.03);
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
            xy_projection_window.draw(point2shape);
            xy_projection_window.draw(pointShape);
        }

        xy_projection_window.draw(axes);
        xy_projection_window.draw(&centerDrawing, 50, sf::Points);
        xy_projection_window.draw(circleDrawing);
        xy_projection_window.draw(circleDrawingUnshuffled);

        xy_projection_window.display();
    }
    sf::RenderWindow xz_projection_window(sf::VideoMode(1000, 1000), "15");
    xz_projection_window.setFramerateLimit(60);
    view.setCenter(0, 0);
    view.setSize(40, 40);
    sf::CircleShape circleXZDrawing(radius);
    circleDrawing.setScale((1 / alpha[0]), 1 / alpha[2]);

    while (xz_projection_window.isOpen()) {
        sf::Event event;
        xz_projection_window.setView(view);

        while (xz_projection_window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                xz_projection_window.close();
            }
        }

        xz_projection_window.clear(sf::Color::Black);

        sf::Vertex centerDrawing(sf::Vector2f(center.x, center.z), sf::Color::Red);

        circleXZDrawing.setPosition((center.x - radius ), (center.z - radius));
        circleXZDrawing.setFillColor(sf::Color::Transparent);
        circleXZDrawing.setOutlineThickness(0.03);
        circleXZDrawing.setOutlineColor(sf::Color(250, 150, 100));

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
            point2shape.setPosition(points[i].x, points[i].z);
            xz_projection_window.draw(point2shape);
            xz_projection_window.draw(pointShape);
        }

        xz_projection_window.draw(axes);
        xz_projection_window.draw(&centerDrawing, 50, sf::Points);
        xz_projection_window.draw(circleDrawing);
        xz_projection_window.draw(circleDrawingUnshuffled);

        xz_projection_window.display();
    }





    return 0;
}

