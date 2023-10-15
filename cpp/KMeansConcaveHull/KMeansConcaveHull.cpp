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
#include <bits/stdc++.h>

#include "KMeansConcaveHull.hpp"

namespace Clustering
{

bool lineIntersects(const vector<Point>& points);
uint32_t NumTrueBools(const vector<bool>& boolVector);

KmeansConcaveHull::KmeansConcaveHull(const vector<double>& lat, const vector<double>& lon)
    : _mask(lat.size(), true)
{
    _data_set.reserve(lat.size());
    for(int i = 0; i < lat.size(); i++)
    {
        _data_set.push_back(Point(lon[i], lat[i]));
    }

    // Sort dataset by Longitude (x value).
    sort(_data_set.begin(), _data_set.end(),
            [](const Point& left, const Point& right) -> bool {
                return left.x < right.x;
            });
}

KmeansConcaveHull::KmeansConcaveHull(const vector<Point>& dataset)
    : _data_set(dataset)
    , _mask(dataset.size(), true)
{
}

vector<Point> KmeansConcaveHull::calculate(size_t k)
{
    return calculate(_data_set, k);
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
    
    vector<Point> hull; 
    hull.push_back(_points[first_point]);
    vector<Point> test_hull = hull;

    _mask[first_point] = false;

    double prev_angle = 270.0l;
    uint32_t step = 2;
    uint32_t stop = step + k_check;
    
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

            Point next_point = _points[knn[candidate]];

            test_hull = hull;
            test_hull.push_back(next_point);

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
        hull.push_back(hull[0]);
        return hull;
    }
    else
    {
        return recurseCalculate(_points);
    }

    return hull;
}

uint32_t KmeansConcaveHull::getLowestLatitudeIndex(const vector<Point>& _points)
{
    vector<double> temp_lats = getLats(_points);
    vector<double>::iterator it = min_element(begin(temp_lats), end(temp_lats));
    uint32_t index = distance(begin(temp_lats), it);
    return index;
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

vector<uint32_t> KmeansConcaveHull::getKNearest(uint32_t currentPointIndex, size_t k)
{
    vector<uint32_t> base_indices = getMaskedIndices(range(_mask.size()), _mask);

    // TODO: work out why warnings exist.
    vector<Point> masked_data_set = arraySubset<Point>(_data_set, base_indices);

    vector<double> distances = calculateDistances(_data_set[currentPointIndex], masked_data_set);

    // Sort the distances array in non-decending order
    vector<uint32_t> sorted_indices = argsort(distances);

    // Get the index of the lowest K points
    size_t k_check = min(k, sorted_indices.size());
    sorted_indices.resize(k_check);

    // Set the index of the lowest K points to True, rest are false
    vector<uint32_t> kNearest = arraySubset<uint32_t>(base_indices, sorted_indices);
    
    return kNearest;
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

vector<double> KmeansConcaveHull::calculateDistances(Point currentPoint, const vector<Point>& kNearestPoints)
{
    vector<double> distances;
    distances.reserve(kNearestPoints.size());
    for(int i = 0; i < kNearestPoints.size(); i++)
    {
        distances.push_back(Clustering::haversineDistance(currentPoint, kNearestPoints[i]));
    }

    return distances;
}

// TODO: Make this dynamic searchPointIndices and return so we can use it with single input/output.
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
        double heading = Clustering::calculateHeading(currentPoint, searchPoints[i], ref_heading);
        headings.push_back(heading);
    }

    return headings;
}

bool lineIntersects(const vector<Point>& points)
{
    for(int i = 0; i < points.size() - 1; i++)
    {
        Point lineOneFirstPoint(points[i].x, points[i].y);
        Point lineOneSecondPoint(points[i + 1].x, points[i + 1].y);
        for(int j = 0; j < points.size() - 1; j++)
        {
            if(i == j) // Don't check intersection of point with itself
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

}