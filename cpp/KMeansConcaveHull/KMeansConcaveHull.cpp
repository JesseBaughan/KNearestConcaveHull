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

#include <iostream>
#include <stdio.h>
 
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <assert.h>
#include <bits/stdc++.h>

#include "KMeansConcaveHull.hpp"
#include "../LineIntersect.hpp"
#include "../PointInPolygonCheck.h"

namespace Clustering
{

static inline double radians(double degrees)
{
    return (degrees * M_PI)  / 180.0f;
}

static inline double degrees(double radians)
{
    return radians * (180.0 / M_PI);
}

KmeansConcaveHull::KmeansConcaveHull(const vector<double>& lat, const vector<double>& lon)
    : _mask(lat.size(), true)
{
    _data_set.reserve(lat.size());
    for(int i = 0; i < lat.size(); i++)
    {
        _data_set.push_back(Point(lon[i], lat[i]));
    }
    sort(_data_set.begin(), _data_set.end(),
            [](const Point& left, const Point& right) -> bool {
                // sort indices according to corresponding array element
                return left.x < right.x;
            });
}

KmeansConcaveHull::KmeansConcaveHull(const vector<Point>& dataset)
    : _data_set(dataset)
    , _mask(dataset.size(), true)
{
}
/*
vector<vector<double>> KmeansConcaveHull::iterativeCalculate()
{
    //Check our data
    //TODO: do we want to do data init that we do in constuctor here and 
    // therefore take the dataset as input to this function?
    // We could reduce constructing this class over and over if we 
    // can just pass in data. Also having NO class globals would be preferred.
    if (_data_set.size() <= 3)
    {
        cout << "Skipped hull calc <= 3 points" << endl;
        //return something?
    }

    vector<vector<double>> hull;
    hull.reserve(_data_set.size());

    for (auto& k: _prime_k)
    {
        uint32_t k_check = min(static_cast<size_t>(k), _data_set.size());

        uint32_t first_point = getLowestLatitudeIndex();
        uint32_t current_point = first_point;
        
        // This needs fixing and working out how to do in c++
        hull = np.reshape(np.array(self.data_set[first_point, :]), (1, 2))
        test_hull = hull;

        _mask[first_point] = false;

        double prev_angle = 270.0l;
        uint32_t step = 2;
        uint32_t stop = step + k_check;
        bool currentKValueFail = false;
        
        // Not sure what the self mask [self mask] is doing.
        //while ((current_point != first_point) or (step == 2)) && (self._mask[self._mask]) > 0:
        // This is the main logic that is stepping through our dataset and attempting to draw 
        // a concanve hull. This could probably be abstracted into it's own function.
        while (((current_point != first_point) or (step == 2)) && (self._mask[self._mask]) > 0)
        {
            if (step == stop)
            {
                _mask[first_point] = true;
            }
            
            vector<uint32_t> knn = getKNearest(current_point, k_check);

            vector<double> angles = calculateHeadings(current_point, knn, prev_angle);

            vector<uint32_t> candidates = np.argsort(-angles);

            uint32_t i = 0;
            bool invalid_hull = true;
            uint32_t candidate = candidates[i++];

            while (invalid_hull && (i < candidates.size()))
            {
                uint32_t candidate = candidates[i];

                // All of this stuff needs working out.
                uint32_t next_point = np.reshape(self.data_set[knn[candidate]], (1,2));
                test_hull = np.append(hull, next_point, axis=0);

                // How the heck are we going to work out if there is a collision of two lines??
                line = LineString(test_hull);
                bool invalid_hull = !line.is_simple;
                i += 1;
            }

            if (invalid_hull)
            {
                currentKValueFail = true;
                break;
            }

            vector<double> prev_angle = calculateHeadings(knn[candidate], np.array([current_point]));
            uint32_t current_point = knn[candidate];
            // Work out how to store this data.
            hull = test_hull;
            
            _mask[current_point] = false;
            step += 1;
        }

        if(currentKValueFail)
        {
            continue;
        }
        
        uint32_t count = 0;
        size_t total = _data_set.size();

        // Check if all our our data set is within our candidate hull 
        for (int ix = 0; ix < total; ix++)
        {
            // TODO: We are going to have to write a contained check.
            if (__contained_check(hull, _data_set[ix, :]))
            {
                count += 1;
            }
            else 
            {
                break;
            }
        }
        
        // If all of our points are within our hull then we can return the hull (success)
        // Else we need to go and try the process again with a new K value.
        if (count == total)
        {
            hull = np.append(hull, [hull[0]], axis=0);
            break; //we have found a solution.
        }
        else
        {
            continue;
        }
    }

    return hull;
}
*/

bool lineIntersects(const vector<Point>& points)
{
    for(int i = 0; i < points.size() - 1; i++)
    {
        Point lineOneFirstPoint(points[i].x, points[i].y);
        Point lineOneSecondPoint(points[i + 1].x, points[i + 1].y);
        for(int j = 0; j < points.size() - 1; j++)
        {
            if(i == j) // don't check intersection of point with itself
            {
                continue;
            }
            else
            {
                Point lineTwoFirstPoint(points[j].x, points[j].y);
                Point lineTwoSecondPoint(points[j + 1].x, points[j + 1].y);
                bool intersects = lineLineIntersection(lineOneFirstPoint, lineOneSecondPoint, 
                                                            lineTwoFirstPoint, lineTwoSecondPoint);
                if(intersects)
                {
                    return true;
                }
            }
        }
    }

    return false;
}

uint32_t NumTrueBools(const vector<bool>& boolVector)
{
    uint32_t numTrueBools = 0;

    // TODO: USE STL COUNT_IF ALGORIGHTM?
    for(const auto& value : boolVector)
    {
        if(value == true)
        {
            numTrueBools++;
        }
    }

    return numTrueBools;
}

vector<Point> KmeansConcaveHull::calculate(const vector<Point>& _points, size_t k)
{
    if (_points.size() <= 3)
    {
        cout << "Skipped hull calc <= 3 points" << endl;
        //return something?
    }
    
    uint32_t k_check = min(k, _points.size());

    uint32_t first_point = getLowestLatitudeIndex(_points);
    uint32_t current_point = first_point;
    
    // This needs fixing and working out how to do in c++
    //hull = np.reshape(np.array(self.data_set[first_point, :]), (1, 2))
    //test_hull = hull;
    vector<Point> hull; 
    hull.push_back(_points[first_point]);
    vector<Point> test_hull = hull;

    _mask[first_point] = false;

    double prev_angle = 270.0l;
    uint32_t step = 2;
    uint32_t stop = step + k_check;
    
    //while ((current_point != first_point) or (step == 2)) && (self._mask[self._mask]) > 0:
    while ((((current_point != first_point) || (step == 2)) && NumTrueBools(_mask)) > 0)
    {
        if (step == stop)
        {
            _mask[first_point] = true;
        }
        
        vector<uint32_t> knn = getKNearest(current_point, k_check);

        vector<double> angles = calculateHeadings(current_point, knn, prev_angle);

        //TODO: why do we negate the array here?
        NegateArray<double>(angles);
        vector<uint32_t> candidates = argsort<double>(angles);

        uint32_t i = 0;
        bool invalid_hull = true;
        uint32_t candidate = 0;

        while (invalid_hull && (i < candidates.size()))
        {
            candidate = candidates[i];

            // All of this stuff needs working out.
            //uint32_t next_point = np.reshape(self.data_set[knn[candidate]], (1,2));
            Point next_point = _points[knn[candidate]];

            test_hull = hull;
            test_hull.push_back(next_point);

            // How the heck are we going to work out if there is a collision of two lines??
            invalid_hull = lineIntersects(test_hull);
            i += 1;
        }

        if (invalid_hull)
        {
            return recurseCalculate(_points);
        }

        vector<uint32_t> current_point_arr = {current_point};
        vector<double> prev_angle_arr = calculateHeadings(knn[candidate], current_point_arr);
        prev_angle = prev_angle_arr[0];
        current_point = knn[candidate];
        // Work out how to store this data.
        hull = test_hull;
        
        _mask[current_point] = false;
        step += 1;
    }
    
    uint32_t count = 0;
    size_t total = _points.size();

    for (int index = 0; index < total; index++)
    {
        Point point{_points[index].x, _points[index].y};
        if (isInside(point, hull))
        {
            count += 1;
        }
        else 
        {
            break;
        }
    }
    
    if (count == total)
    {
        //TODO: confirm this operation matches python.
        hull.push_back(hull[0]);
        return hull;
    }
    else
    {
        return recurseCalculate(_points);
    }

    vector<Point> poo(5);
    return poo;
}

int KmeansConcaveHull::getNextK()
{
    if (_current_prime_index++ < _prime_k.size())
    {
        return _prime_k[_current_prime_index];
    }
    else 
    {
        return -1;
    }
}

vector<Point> KmeansConcaveHull::recurseCalculate(const vector<Point>& points, uint32_t k)
{
    int next_k = getNextK();
    if (next_k == -1)
    {   
        vector<Point> empty;
        return empty;
    }

    Clustering::KmeansConcaveHull hullCalc(points);

    return hullCalc.calculate(points, next_k);
}

double KmeansConcaveHull::haversineDistance(const Point first, const Point second)
{
    const double earths_radius = 6371000.0f;

    // Get the difference between our two points then radians the difference into radians
    const double lat_delta = radians(second.y - first.y);
    const double lon_delta = radians(second.x - first.x);

    const double converted_lat1 = radians(first.y);
    const double converted_lat2 = radians(second.y);

    const double a =
        pow(sin(lat_delta / 2), 2) + cos(converted_lat1) * cos(converted_lat2) * pow(sin(lon_delta / 2), 2);

    const double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    const double d = earths_radius * c;

    return d;
}

vector<double> KmeansConcaveHull::calculateDistances(Point currentPoint, const vector<Point>& kNearestPoints)
{
    vector<double> distances;
    distances.reserve(kNearestPoints.size());
    for(int i = 0; i < kNearestPoints.size(); i++)
    {
        distances.push_back(haversineDistance(currentPoint, kNearestPoints[i]));
    }

    return distances;
}

vector<double> KmeansConcaveHull::getLats(const vector<Point>& coords)
{
    vector<double> lats;
    lats.reserve(coords.size());
    for(auto& coord: coords)
    {
        lats.push_back(coord.y);
    }

    return lats;
}

uint32_t KmeansConcaveHull::getLowestLatitudeIndex(const vector<Point>& _points)
{
    vector<double> temp_lats = getLats(_points);
    vector<double>::iterator it = min_element(begin(temp_lats), end(temp_lats));
    uint32_t index = distance(begin(temp_lats), it);
    return index;
}

vector<uint32_t> KmeansConcaveHull::getMaskedIndices(const vector<uint32_t>& input_array, const vector<bool>& mask)
{
    vector<uint32_t> masked_array;
    masked_array.reserve(input_array.size());
    for(int i = 0; i < input_array.size(); i++)
    {
        if(mask[i] == true)
        {
            masked_array.push_back(input_array[i]);
        }
    }

    return masked_array;
}

// Returns a range from 0 to size
vector<uint32_t> range(size_t size)
{
    vector<uint32_t> output_array(size, 0);
    iota(output_array.begin(), output_array.end(), 0);
    return output_array;
}

vector<uint32_t> KmeansConcaveHull::getKNearest(uint32_t currentPointIndex, size_t k)
{
    // Harcoded for testing purposes.
    //_mask[1] = false;
    vector<uint32_t> base_indices = getMaskedIndices(range(_mask.size()), _mask);

    vector<Point> masked_data_set = arraySubset<Point>(_data_set, base_indices);

    vector<double> distances = calculateDistances(_data_set[currentPointIndex], masked_data_set);

    //Sort the distances array in non-decending order
    vector<uint32_t> sorted_indices = argsort(distances);

    //Get the index of the lowest K points
    size_t k_check = min(k, sorted_indices.size());
    sorted_indices.resize(k_check);

    //Set the index of the lowest K points to True, rest are false
    vector<uint32_t> kNearest = arraySubset<uint32_t>(base_indices, sorted_indices);
    
    return kNearest;
}

/**
 * Argsort(currently support ascending sort)
 * @tparam T array element type
 * @param array input array
 * @return indices w.r.t sorted array
 */
template<typename T>
vector<uint32_t> KmeansConcaveHull::argsort(const vector<T> &array) {
    vector<uint32_t> indices(array.size());
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(),
            [&array](int left, int right) -> bool {
                // sort indices according to corresponding array element
                return array[left] < array[right];
            });

    return indices;
}

double KmeansConcaveHull::calculateHeading(Point reference, Point target, double ref_heading)
{
    if ((ref_heading < 0.0f) || (ref_heading >= 360.0f))
    {
        assert(false);
    }
    
    double referencePointLat_rads = radians(reference.y);
    double referencePointLon_rads = radians(reference.x);
    double targetPointLat_rads = radians(target.y);
    double targetPointLon_rads = radians(target.x);

    double lon_dif = targetPointLon_rads - referencePointLon_rads;

    float y = sin(lon_dif) * cos(targetPointLat_rads);

    float x = cos(referencePointLat_rads) * sin(targetPointLat_rads) - sin(referencePointLat_rads) *
            cos(targetPointLat_rads) * cos(lon_dif);

    float bearing = fmod((degrees(atan2(y, x)) + 360.0l), 360.0l) - ref_heading;

    if(bearing < 0.0l)
    {
        bearing += 360.0l;
    }
    return bearing;
}

//Make this dynamic searchPointIndices and return so we can use it with single input/output
vector<double> KmeansConcaveHull::calculateHeadings(uint32_t currentPointIndex, 
                                                        const vector<uint32_t>& searchPointsIndicies, 
                                                        double ref_heading)
{
    // Get the subset of points we are calculating heading to 
    vector<Point> searchPoints = arraySubset(_data_set, searchPointsIndicies);
    Point currentPoint = _data_set[currentPointIndex];

    vector<double> headings;
    headings.reserve(searchPoints.size());
    for(int i = 0; i < searchPoints.size(); i++)
    {
        double heading = calculateHeading(currentPoint, searchPoints[i], ref_heading);
        headings.push_back(heading);
    }

    return headings;
}

bool containedCheck(const vector<Point>& hull, Point point)
{
    return true;
}

}