/**
 * @file KMeansConcaveHull.cpp
 *
 * @brief A description of the module’s purpose.
 *
 * @author Name - email address
 *
 * @par 
 *
 */
 
#include <math.h>
#include <cmath>

#include "KMeansConcaveHull.hpp"

// convert our passed value to radians_t
inline float convert(const float angle) 
{
    return angle * (M_PI / 180);
}

float KmeansConcaveHull::havesineDistance(lat_lon_coord first, lat_lon_coord second)
{
    const float earths_radius = 6371;

    // Get the difference between our two points then convert the difference into radians
    const float lat_delta = convert(second.Lat - first.Lat);
    const float lon_delta = convert(second.Lon - first.Lon);

    const float converted_lat1 = convert(first.Lat);
    const float converted_lat2 = convert(second.Lat);

    const float a =
        pow(sin(lat_delta / 2), 2) + cos(converted_lat1) * cos(converted_lat2) * pow(sin(lon_delta / 2), 2);

    const float c = 2 * atan2(sqrt(a), sqrt(1 - a));
    const float d = earths_radius * c;

    return d;
}

