#pragma once

struct lat_lon_coord
{
    double Lat;
    double Lon;

    lat_lon_coord(double lat, double lon)
        : Lat(lat)
        , Lon(lon) {}
    
    lat_lon_coord()
    : Lat(0)
    , Lon(0) {}
};

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