/**
 * @file GeometryUtils.h
 *
 * @brief Geometry utility functions used to determine the concave hull.
 * 
 * @author Jesse B - email address
 *
 */

#ifndef GEOMETRY_UTILS_H_
#define GEOMETRY_UTILS_H_

#include <vector>

namespace Clustering
{
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

/**
 * @brief Can be used to check the inersection of two lines
 * defined by their start/end points.
* */
bool lineLineIntersection(Point A, Point B, Point C, Point D);

bool isInside(const Point &point, const std::vector<Point> &points_list);

bool PointLiesOnLine(Point intersectPoint, Point A, Point B);

}

#endif /* GEOMETRY_UTILS_H_ */

/*** end of file ***/