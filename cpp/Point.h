#pragma once

// This pair is used to store the X and Y
// coordinates of a point respectively
struct Point 
{
    double x;
    double y;

    Point()
        : x(0.0l)
        , y(0.0l) {}

    Point(double _x, double _y)
        : x(_x)
        , y(_y) {}
};