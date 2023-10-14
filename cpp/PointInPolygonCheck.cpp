#include <limits>
#include "PointInPolygonCheck.h"

int Is_Left(const Point &p0, const Point &p1, const Point &point)
{
    return ((p1.x - p0.x) * (point.y - p0.y) -
            (point.x - p0.x) * (p1.y - p0.y));
}

template<typename T>
int Is_Equal(const T val1, const T val2)
{
   double diff =  val1 - val2;
   return (diff < std::numeric_limits<double>::epsilon() ) && (-diff < std::numeric_limits<double>::epsilon());
}

int Is_Intersecting(const Point &p1, const Point& p2)
{
    bool xIntersects = Is_Equal(p1.x, p2.x);
    bool yIntersects = Is_Equal(p1.y, p2.y);
    return xIntersects || yIntersects;
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
        bool p1Intersecting = Is_Intersecting(point, point1);
        bool p2Intersecting = Is_Intersecting(point, point2);
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
                    winding_number += 1.0f;                               // Have a valid up intersect
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
                    winding_number -= 1.0f;                               // Have a valid down intersect
                }
            }
        }
    }

    return (winding_number > 0.25f) || (winding_number < -0.25);
}