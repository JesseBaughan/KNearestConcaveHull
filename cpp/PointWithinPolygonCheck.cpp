#include "PointWithinPolygonCheck.h"

//TODO: generalise to just take x/y not Lat/Lon
//X = Lat, Y = Lon at the current stage.
int pointLiesWithinPolygon(vector<lat_lon_coord>& points, float testx, float testy)
{
    int i, j, c = 0;
    for (i = 0, j = (points.size() - 1); i < points.size(); j = i++) {
        if ( ((points[i].Lon > testy) != (points[j].Lon > testy)) &&
        (testx < (points[j].Lat - points[i].Lat) * (testy - points[i].Lon) / (points[j].Lon - points[i].Lon) + points[i].Lat) )
        c = !c;
    }
    return c;
}