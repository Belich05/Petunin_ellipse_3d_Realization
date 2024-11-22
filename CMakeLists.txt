cmake_minimum_required(VERSION 3.26)
project(kolok2)
include (FetchContent)
set(BUILD_SHARED_LIBS OFF)
FETCHCONTENT_DECLARE (SFML GIT_REPOSITORY https://github.com/SFML/SFML.git GIT_TAG 2.5.1)
FETCHCONTENT_MAKEAVAILABLE (SFML)
set (CMAKE_CXX_STANDARD 17)
add_executable(kolok2
        main.cpp
        Point.h
        PlaneEquation.h
        Sphere.h
        Circle.h
        LineEquation.h
        Rectangle.h
        Cube.h
        MathTools.h
)
target_link_libraries (kolok2 sfml-graphics)