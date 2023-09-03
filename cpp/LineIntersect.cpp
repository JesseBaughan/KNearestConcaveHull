#include <vector>
#include <cmath>

#include "LineIntersect.hpp"

#define FLT_MAX 999999

/*
bool LineIntersects(const std::vector<lat_lon_coord>& points)
{
    for(int i = 0; i < points.size() - 1; i++)
    {
        Point lineOneFirstPoint(points[i].Lat, points[i].Lon);
        Point lineOneSecondPoint(points[i + 1].Lat, points[i + 1].Lon);
        for(int j = 0; j < points.size() - 1; j++)
        {
            if(i == j)
            {
                continue;
            }
            else
            {
                Point lineTwoFirstPoint(points[j].Lat, points[j].Lon);
                Point lineTwoSecondPoint(points[j + 1].Lat, points[j + 1].Lon);
                Point intersectPoint = lineLineIntersection(lineOneFirstPoint, lineOneSecondPoint, 
                                                            lineTwoFirstPoint, lineTwoSecondPoint);
            }
        }
    }
}

*/

bool PointLiesOnLine(Point intersectPoint, Point A, Point B)
{
    bool xPointOnLine = (intersectPoint.first >= fmin(A.first, B.first)) && 
                        (intersectPoint.first <= fmax(A.first, B.first));

    bool yPointOnLine = (intersectPoint.second >= fmin(A.second, B.second)) && 
                        (intersectPoint.second <= fmax(A.second, B.second));

    return xPointOnLine && yPointOnLine;
}

bool LineLineIntersection(Point A, Point B, Point C, Point D)
{
    // Line AB represented as a1x + b1y = c1
    double a1 = B.second - A.second;
    double b1 = A.first - B.first;
    double c1 = (a1 * A.first) + (b1 * A.second);

    // Line CD represented as a2x + b2y = c2
    double a2 = D.second - C.second;
    double b2 = C.first - D.first;
    double c2 = (a2 * C.first) + (b2 * C.second);

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
        return (PointLiesOnLine(Point(x, y), A, B) && PointLiesOnLine(Point(x, y), C, D));
    }
}