
/**
 * @file LineIntersect.cpp
 *
 * @brief Can be used to check the inersection of two lines
 * defined by their start/end points.
 *
 * @author Jesse B - email address
 *
 */

// This pair is used to store the X and Y
// coordinates of a point respectively
struct Point 
{
    double first;
    double second;

    Point(double _first, double _second)
        : first(_first)
        , second(_second) {}
};

bool lineLineIntersection(Point A, Point B, Point C, Point D);