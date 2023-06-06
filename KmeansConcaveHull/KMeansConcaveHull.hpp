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
#include <array>

namespace Clustering
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
        KmeansConcaveHull(const std::vector<float>& lat, const std::vector<float>& lon)
            : _lat(lat)
            , _lon(lon)
            , _indices(_lat.size(), true)
        {
        }

        ~KmeansConcaveHull() {};

        KmeansConcaveHull(const KmeansConcaveHull&) = delete;
        KmeansConcaveHull& operator = (const KmeansConcaveHull&) = delete;

        std::vector<std::vector<float>> calculate(const std::vector<std::vector<float>>& points, uint32_t k = 0);
        uint32_t getLowestLatitudeIndex();

        std::vector<bool> get_indices()    const { return _indices; }

    private:
        const std::array<float, 18> prime_k = {3,  7, 13, 19, 29, 37,
                                               43, 53, 61, 71, 79, 89,
                                               97, 101, 107, 113, 131, 139};
        std::vector<lat_lon_coord> _data_set;
        std::vector<float> _lat;
        std::vector<float> _lon;
        std::vector<bool> _indices;

        uint32_t prime_ix{0};

        std::vector<float> calculateDistances(lat_lon_coord currentPoint, const std::vector<lat_lon_coord>& kNearestPoints);
        float haversineDistance(lat_lon_coord first, lat_lon_coord second);

        std::vector<bool> getKNearest(uint32_t currentPointIndex, uint32_t k);

        float getNextK();

        std::vector<float> calculateHeadings(lat_lon_coord currentPointIndex, 
                                             const std::vector<lat_lon_coord>& searchPoints, 
                                             float ref_heading=0.0f);

        bool containedCheck(const std::vector<lat_lon_coord>& hull, lat_lon_coord point);

        // TODO: refacto to NOT use recursion - it is not efficient.
        void recurseCalculate(const std::vector<lat_lon_coord>& points, uint32_t k = 3);

        template<typename Type>
        std::vector<Type> KmeansConcaveHull::getMaskedArray(const std::vector<Type>& input_array, std::vector<bool>& mask);

    };
}

#endif /* KMEANS_CONCAVE_HULL_H */

/*** end of file ***/