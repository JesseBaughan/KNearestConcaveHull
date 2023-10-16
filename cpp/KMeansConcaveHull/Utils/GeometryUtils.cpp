/**
 * @file GeometryUtils.cpp
 *
 * @brief Geometry utility functions used to determine the concave hull.
 * 
 * @author Jesse B - email address
 *
 */

#include <limits>
#include <vector>
#include <cmath>
#include <assert.h>

#include "GeometryUtils.h"

#define FLT_MAX 999999

using namespace std;

namespace Clustering
{

int Is_Left(const Point &p0, const Point &p1, const Point &point);

int Is_PointEqual(const Point &p1, const Point& p2);

bool lineLineIntersection(Point A, Point B, Point C, Point D)
{
    // Line AB represented as a1x + b1y = c1
    double a1 = B.y - A.y;
    double b1 = A.x - B.x;
    double c1 = (a1 * A.x) + (b1 * A.y);

    // Line CD represented as a2x + b2y = c2
    double a2 = D.y - C.y;
    double b2 = C.x - D.x;
    double c2 = (a2 * C.x) + (b2 * C.y);

    double determinant = (a1 * b2) - (a2 * b1);

    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of FLT_MAX
        return false;
    }
    else
    {
        double x = (b2 * c1 - b1 * c2) / determinant;
        double y = (a1 * c2 - a2 * c1) / determinant;
        return (pointLiesOnLine(Point(x, y), A, B) && pointLiesOnLine(Point(x, y), C, D));
    }
}

bool pointLiesOnLine(Point intersectPoint, Point A, Point B)
{
    bool xPointOnLine = (intersectPoint.x > fmin(A.x, B.x)) && 
                        (intersectPoint.x < fmax(A.x, B.x));
    bool yPointOnLine = (intersectPoint.y > fmin(A.y, B.y)) && 
                        (intersectPoint.y < fmax(A.y, B.y));

    return xPointOnLine && yPointOnLine;
}

// Paper explaining this algorithm can be found here: 
// https://www.engr.colostate.edu/~dga/documents/papers/point_in_polygon.pdf
// Source code was pulled from this page: 
// https://forums.codeguru.com/showthread.php?497679-To-check-if-a-point-is-inside-a-polygon 
bool isInside(const Point &point, const std::vector<Point> &points_list)
{
    // The winding number counter.
    float winding_number = 0.0f;

    // Loop through all edges of the polygon.
    typedef std::vector<Point>::size_type size_type;

    size_type size = points_list.size();

    for (size_type i = 0; i < size; ++i)             // Edge from point1 to points_list[i+1]
    {
        Point point1(points_list[i]);
        Point point2;

        // Wrap?
        if (i == (size - 1))
        {
            point2 = points_list[0];
        }
        else
        {
            point2 = points_list[i + 1];
        }

        // We consider a point that is the same as a vertice of the polygon 
        // as inside the polygon.
        bool p1Intersecting = Is_PointEqual(point, point1);
        bool p2Intersecting = Is_PointEqual(point, point2);
        if(p1Intersecting || p2Intersecting)
        {
            return true;
        }

        if (point1.y <= point.y)                                    // start y <= point.y
        {
            if (point2.y > point.y)                                 // An upward crossing
            {
                if (Is_Left(point1, point2, point) > 0)             // Point left of edge
                {
                    winding_number += 1.0f;                         // Have a valid up intersect
                }
            }
        }
        else
        {
            // start y > point.y (no test needed)
            if (point2.y <= point.y)                                // A downward crossing
            {
                if (Is_Left(point1, point2, point) < 0)             // Point right of edge
                {
                    winding_number -= 1.0f;                         // Have a valid down intersect
                }
            }
        }
    }

    return (winding_number > 0.25f) || (winding_number < -0.25);
}

int Is_Left(const Point &p0, const Point &p1, const Point &point)
{
    return ((p1.x - p0.x) * (point.y - p0.y) -
            (point.x - p0.x) * (p1.y - p0.y));
}

int Is_PointEqual(const Point &p1, const Point& p2)
{
    bool xIntersects = isEqual(p1.x, p2.x);
    bool yIntersects = isEqual(p1.y, p2.y);
    return xIntersects || yIntersects;
}

int isEqual(const double val1, const double val2)
{
   double diff =  val1 - val2;
   return (diff < std::numeric_limits<double>::epsilon() ) && (-diff < std::numeric_limits<double>::epsilon());
}

int isEqual(const float val1, const float val2)
{
   double diff =  val1 - val2;
   return (diff < std::numeric_limits<float>::epsilon() ) && (-diff < std::numeric_limits<float>::epsilon());
}

double haversineDistance(const Point first, const Point second)
{
    const double earths_radius = 6371000.0f;

    // Get the difference between our two points then radians the difference into radians
    const double lat_delta = radians(second.y - first.y);
    const double lon_delta = radians(second.x - first.x);

    const double converted_lat1 = radians(first.y);
    const double converted_lat2 = radians(second.y);

    const double a =
        pow(sin(lat_delta / 2), 2) + cos(converted_lat1) * cos(converted_lat2) * pow(sin(lon_delta / 2), 2);

    const double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    const double d = earths_radius * c;

    return d;
}

double calculateHeading(Point reference, Point target, double ref_heading)
{
    if ((ref_heading < 0.0f) || (ref_heading >= 360.0f))
    {
        assert(false);
    }
    
    double referencePointLat_rads = radians(reference.y);
    double referencePointLon_rads = radians(reference.x);
    double targetPointLat_rads = radians(target.y);
    double targetPointLon_rads = radians(target.x);

    double lon_dif = targetPointLon_rads - referencePointLon_rads;

    float y = sin(lon_dif) * cos(targetPointLat_rads);

    float x = cos(referencePointLat_rads) * sin(targetPointLat_rads) - sin(referencePointLat_rads) *
            cos(targetPointLat_rads) * cos(lon_dif);

    float bearing = fmod((degrees(atan2(y, x)) + 360.0l), 360.0l) - ref_heading;

    if(bearing < 0.0l)
    {
        bearing += 360.0l;
    }
    return bearing;
}

}

/*** end of file ***/