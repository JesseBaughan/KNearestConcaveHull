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

    class KmeansConcaveHull
    {
    public:
        KmeansConcaveHull(const vector<double>& lat, const vector<double>& lon);
        explicit KmeansConcaveHull(const vector<Point>& dataset);
        ~KmeansConcaveHull() {};

        KmeansConcaveHull(const KmeansConcaveHull&) = delete;
        KmeansConcaveHull& operator = (const KmeansConcaveHull&) = delete;

        vector<Point> calculate(const vector<Point>& points, size_t k);

        // These are temporarily public whilst we perform testing
        vector<bool> get_mask()    const { return _mask; }
        uint32_t getLowestLatitudeIndex(const vector<Point>& points);
        vector<uint32_t> getKNearest(uint32_t currentPointIndex, size_t k = 3);

        vector<double> calculateHeadings(uint32_t currentPointIndex, 
                                             const vector<uint32_t>& searchPoints, 
                                             double ref_heading=0.0l);

        // Public for development purposes
        vector<Point> _data_set;

    private:
        const array<int, 18> _prime_k = {3,  7, 13, 19, 29, 37,
                                         43, 53, 61, 71, 79, 89,
                                         97, 101, 107, 113, 131, 139};
        size_t _current_prime_index{0};

        vector<bool> _mask;

        uint32_t _prime_ix{0};
        size_t _k{3};

        vector<double> calculateDistances(Point currentPoint, 
                                              const vector<Point>& kNearestPoints);

        double haversineDistance(Point first, Point second);

        int getNextK();

        double calculateHeading(Point reference, Point target, double ref_heading);

        vector<Point> recurseCalculate(const vector<Point>& points, uint32_t k = 3);

        vector<uint32_t> getMaskedIndices(const vector<uint32_t>& input_array, 
                                                                  const vector<bool>& mask);
        vector<double> getLats(const vector<Point>& coords);
    };

   template<typename T>
   vector<uint32_t> argsort(const vector<T> &array);

   vector<uint32_t> range(size_t size);

   template<typename T>
   void NegateArray(vector<T>& inputVector)
   {
        for(int i = 0; i < inputVector.size(); i++)
        {
            inputVector[i] = inputVector[i] * -1;
        }
   }

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
}

#endif /* KMEANS_CONCAVE_HULL_H */

/*** end of file ***/