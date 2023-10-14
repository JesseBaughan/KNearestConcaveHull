#include "PointInPolygonCheck.h"

int Is_Left(const Point &p0,
            const Point &p1,
            const Point &point)
{
    return ((p1.x - p0.x) * (point.y - p0.y) -
            (point.x - p0.x) * (p1.y - p0.y));
}

// Paper explaining this algorithm can be found here: 
// https://www.engr.colostate.edu/~dga/documents/papers/point_in_polygon.pdf
// Source code was pulled from this page: 
// https://forums.codeguru.com/showthread.php?497679-To-check-if-a-point-is-inside-a-polygon 
bool isInside(const Point &point, const std::vector<Point> &points_list)
{
    // The winding number counter.
    int winding_number = 0;

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

        if (point1.y <= point.y)                                    // start y <= point.y
        {
            if (point2.y > point.y)                                 // An upward crossing
            {
                if (Is_Left(point1, point2, point) > 0)             // Point left of edge
                {
                    ++winding_number;                               // Have a valid up intersect
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
                    --winding_number;                               // Have a valid down intersect
                }
            }
        }
    }

    return (winding_number != 0);
}