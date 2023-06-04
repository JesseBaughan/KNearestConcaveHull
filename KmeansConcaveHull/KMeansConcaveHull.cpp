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

lat_lon_coord KmeansConcaveHull::getLowestLatitudeIndex(const std::vector<lat_lon_coord> &points)
{
    // Search array of points we are trying to 
    // calcualte a concave hull on for lowest lat point
    // TODO: What search algorithm do we use here?
}

std::vector<lat_lon_coord> KmeansConcaveHull::getKNearest(uint32_t currentPointIndex, uint32_t k)
{
    /*
    ixs = self.indices

    base_indices = np.arange(len(ixs))[ixs]
    distances = self.haversine_distance(self.data_set[ix, :], self.data_set[ixs, :])
    sorted_indices = np.argsort(distances)

    k_check = min(k, len(sorted_indices))
    k_nearest = sorted_indices[range(k_check)]
    return base_indices[k_nearest]
    */
}

float KmeansConcaveHull::getNextK()
{
    if (prime_ix < prime_k.size())
    {
        return prime_k[prime_ix];
    }
    else
    {
        return -1.0f;
    }
}
