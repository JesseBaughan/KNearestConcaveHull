/**
 * @file KMeansConcaveHull.hpp
 *
 * @brief A description of the module’s purpose.
 *
 * @author Name - email address
 *
 * @par 
 *
 */
 
#ifndef KMEANS_CONCAVE_HULL_H
#define KMEANS_CONCAVE_HULL_H

#include <stdint.h>
#include <vector>

namespace
{
    struct lat_lon_coord
    {
        float Lat;
        float Lon;

        lat_lon_coord(float lat, float lon)
            : Lat(lat)
            , Lon(lon) {}
    };

    class KmeansConcaveHull
    {
    public:
        KmeansConcaveHull();

        ~KmeansConcaveHull();

        std::vector<std::vector<float>> calculate(const std::vector<std::vector<float>>& points, uint32_t k = 0);

    private:
        uint32_t getNextK();

        float havesineDistance(lat_lon_coord first, lat_lon_coord second);

        lat_lon_coord getLowestLatitudeIndex(const std::vector<lat_lon_coord>& points);

        std::vector<lat_lon_coord> getKNearest(uint32_t currentPointIndex, uint32_t k);

        std::vector<float> calculateHeadings(lat_lon_coord currentPointIndex, 
                                             const std::vector<lat_lon_coord>& searchPoints, 
                                             float ref_heading=0.0f);

        bool containedCheck(const std::vector<lat_lon_coord>& hull, lat_lon_coord point);

        // TODO: refacto to NOT use recursion - it is not efficient.
        void recurse_calculate(const std::vector<lat_lon_coord>& points, uint32_t k = 3);

    };
}

#endif /* KMEANS_CONCAVE_HULL_H */

/*** end of file ***/