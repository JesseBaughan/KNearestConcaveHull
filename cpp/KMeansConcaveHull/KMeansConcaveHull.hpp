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

#include "../Point.h"

using namespace std;

namespace Clustering
{

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
        explicit KmeansConcaveHull(const vector<lat_lon_coord>& dataset);
        ~KmeansConcaveHull() {};

        KmeansConcaveHull(const KmeansConcaveHull&) = delete;
        KmeansConcaveHull& operator = (const KmeansConcaveHull&) = delete;

        vector<lat_lon_coord> calculate(const vector<lat_lon_coord>& points, size_t k);

        //vector<vector<double>> KmeansConcaveHull::iterativeCalculate();

        // These are temporarily public whilst we perform testing
        vector<bool> get_mask()    const { return _mask; }
        uint32_t getLowestLatitudeIndex(const vector<lat_lon_coord>& points);
        vector<uint32_t> getKNearest(uint32_t currentPointIndex, size_t k = 3);

        vector<double> calculateHeadings(uint32_t currentPointIndex, 
                                             const vector<uint32_t>& searchPoints, 
                                             double ref_heading=0.0l);

        // Public for development purposes
        vector<lat_lon_coord> _data_set;

    private:
        const array<int, 18> _prime_k = {3,  7, 13, 19, 29, 37,
                                         43, 53, 61, 71, 79, 89,
                                         97, 101, 107, 113, 131, 139};
        size_t _current_prime_index{0};

        vector<bool> _mask;

        uint32_t _prime_ix{0};
        size_t _k{3};

        vector<double> calculateDistances(lat_lon_coord currentPoint, 
                                              const vector<lat_lon_coord>& kNearestPoints);

        double haversineDistance(lat_lon_coord first, lat_lon_coord second);

        int getNextK();

        double calculateHeading(lat_lon_coord reference, lat_lon_coord target, double ref_heading);

        bool containedCheck(const vector<lat_lon_coord>& hull, lat_lon_coord point);

        // TODO: refacto to NOT use recursion - it is not efficient.
        vector<lat_lon_coord> recurseCalculate(const vector<lat_lon_coord>& points, uint32_t k = 3);

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

        vector<double> getLats(const vector<lat_lon_coord>& coords);
    };
}

#endif /* KMEANS_CONCAVE_HULL_H */

/*** end of file ***/