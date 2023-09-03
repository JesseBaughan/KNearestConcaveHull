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

using namespace std;

namespace Clustering
{
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

   template<typename T>
   void NegateArray(vector<T>& inputVector)
   {
        for(int i = 0; i < inputVector.size(); i++)
        {
            inputVector[i] = inputVector[i] * -1;
        }
   }

    class KmeansConcaveHull
    {
    public:
        KmeansConcaveHull(const vector<double>& lat, const vector<double>& lon);
        ~KmeansConcaveHull() {};

        KmeansConcaveHull(const KmeansConcaveHull&) = delete;
        KmeansConcaveHull& operator = (const KmeansConcaveHull&) = delete;

        vector<vector<double>> calculate(size_t k);

        //vector<vector<double>> KmeansConcaveHull::iterativeCalculate();

        // These are temporarily public whilst we perform testing
        vector<bool> get_mask()    const { return _mask; }
        uint32_t getLowestLatitudeIndex();
        vector<uint32_t> getKNearest(uint32_t currentPointIndex, size_t k = 3);

        vector<double> calculateHeadings(uint32_t currentPointIndex, 
                                             const vector<uint32_t>& searchPoints, 
                                             double ref_heading=0.0l);

    private:
        const array<double, 18> _prime_k = {3,  7, 13, 19, 29, 37,
                                               43, 53, 61, 71, 79, 89,
                                               97, 101, 107, 113, 131, 139};
        vector<lat_lon_coord> _data_set;
        vector<double> _lat;
        vector<double> _lon;
        vector<bool> _mask;

        uint32_t _prime_ix{0};
        size_t _k{3};

        vector<double> calculateDistances(lat_lon_coord currentPoint, 
                                              const vector<lat_lon_coord>& kNearestPoints);

        double haversineDistance(lat_lon_coord first, lat_lon_coord second);

        double getNextK();

        double calculateHeading(lat_lon_coord reference, lat_lon_coord target, double ref_heading);

        bool containedCheck(const vector<lat_lon_coord>& hull, lat_lon_coord point);

        // TODO: refacto to NOT use recursion - it is not efficient.
        void recurseCalculate(const vector<lat_lon_coord>& points, uint32_t k = 3);

        vector<uint32_t> getMaskedIndices(const vector<uint32_t>& input_array, 
                                                                  const vector<bool>& mask);

        template<typename T>
        vector<uint32_t> argsort(const vector<T> &array);

        template<typename Type>
        vector<Type> arraySubset(const vector<Type>& input_array, const vector<uint32_t>& indexes)
        {
            vector<Type> output_array;
            output_array.reserve(indexes.size());
            for(int i = 0; i < indexes.size(); i++)
            {
                output_array.push_back(input_array[indexes[i]]);
            }

            return output_array;
        }

        vector<double> getLats(vector<lat_lon_coord>& coords);
    };
}

#endif /* KMEANS_CONCAVE_HULL_H */

/*** end of file ***/