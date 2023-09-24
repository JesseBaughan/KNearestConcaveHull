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