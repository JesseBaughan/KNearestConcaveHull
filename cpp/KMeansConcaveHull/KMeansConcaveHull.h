/**
 * @file KMeansConcaveHull.h
 *
 * @brief Attempts to find the concave hull of a set of 
 * lat/lon coordinates. 
 *
 * @author Jesse B 
 */

#ifndef KMEANS_CONCAVE_HULL_H_
#define KMEANS_CONCAVE_HULL_H_

#include <stdint.h>
#include <vector>
#include <array>

#include "Utils/GeometryUtils.h"
#include "Utils/ArrayUtils.h"

using namespace std;

namespace Clustering
{

class KmeansConcaveHull
{
public:
    KmeansConcaveHull(const vector<double>& lat, const vector<double>& lon);
    explicit KmeansConcaveHull(const vector<Point>& dataset);
    ~KmeansConcaveHull() = default;

    vector<Point> calculate(size_t k = 3);

private:
    const array<int, 18> _prime_k = {3,  7, 13, 19, 29, 37,
                                        43, 53, 61, 71, 79, 89,
                                        97, 101, 107, 113, 131, 139};
    size_t _current_prime_index{0};

    vector<Point> _data_set;

    vector<bool> _mask;

    uint32_t _prime_ix{0};
    size_t _k{3};

    vector<Point> calculate(const vector<Point>& points, size_t k);

    vector<bool> getMask()    const { return _mask; }

    uint32_t getLowestLatitudeIndex(const vector<Point>& points);

    vector<uint32_t> getKNearest(uint32_t currentPointIndex, size_t k = 3);

    vector<double> calculateHeadings(uint32_t currentPointIndex, 
                                            const vector<uint32_t>& searchPoints, 
                                            double ref_heading=0.0l);

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

}

#endif /* KMEANS_CONCAVE_HULL_H_ */

/*** end of file ***/