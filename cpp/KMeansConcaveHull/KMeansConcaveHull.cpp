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

#include "KMeansConcaveHull.hpp"

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

    KmeansConcaveHull::KmeansConcaveHull(const std::vector<double>& lat, const std::vector<double>& lon)
        : _lat(lat)
        , _lon(lon)
        , _mask(_lat.size(), true)
    {
        _data_set.reserve(lat.size());
        for(int i = 0; i < lat.size(); i++)
        {
            _data_set.push_back(lat_lon_coord(lat[i], lon[i]));
        }
        std::sort(_data_set.begin(), _data_set.end(),
                [](const lat_lon_coord& left, const lat_lon_coord& right) -> bool {
                    // sort indices according to corresponding array element
                    return left.Lat < right.Lat;
                });
    }

    std::vector<std::vector<double>> KmeansConcaveHull::calculate(const std::vector<std::vector<double>> &points, size_t k)
    {
        /*
        if isinstance(points, np.core.ndarray):
            self.data_set = points
        else:
            raise ValueError("Please provide a [N,2] numpy array")
        self.data_set = np.unique(self.data_set, axis=0)
        
        self.indices = np.ones(self.data_set.shape[0], dtype=bool)

        try:
            if self.data_set.shape[0] <= 3:
                print("Skipped hull calc <= 3 points")
                #timing_logger.info("Skipped hull calc <= 3 points")
                return self.data_set
            
            k_check = min(k, self.data_set.shape[0]) 

            first_point = self.get_lowest_latitude_index(self.data_set)
            current_point = first_point
            
            hull = np.reshape(np.array(self.data_set[first_point, :]), (1, 2))
            test_hull = hull

            self.indices[first_point] = False

            prev_angle = 270
            step = 2
            stop = step + k_check
            
            while ((current_point != first_point) or (step == 2)) and len(self.indices[self.indices]) > 0:
                if step == stop:
                    self.indices[first_point] = True
                
                knn = self.get_k_nearest(current_point, k_check)

                angles = self.calculate_headings(current_point, knn, prev_angle)

                candidates = np.argsort(-angles)

                i = 0
                invalid_hull = True
                
                while invalid_hull and i < len(candidates):
                    candidate = candidates[i]

                    next_point = np.reshape(self.data_set[knn[candidate]], (1,2))
                    test_hull = np.append(hull, next_point, axis=0)
                    line = LineString(test_hull)
                    invalid_hull = not line.is_simple
                    i += 1

                if invalid_hull:
                    return self.recurse_calculate()

                prev_angle = self.calculate_headings(knn[candidate], np.array([current_point]))
                current_point = knn[candidate]
                hull = test_hull
                
                self.indices[current_point] = False
                step += 1
            
            count = 0
            total = self.data_set.shape[0]
            for ix in range(total):
                if self.__contained_check(hull, self.data_set[ix, :]):
                    count += 1
                else:
                    break
            
            if count == total:
                hull = np.append(hull, [hull[0]], axis=0)
                return hull
            else: 
                return self.recurse_calculate()
        except Exception as e:
            print("HullCalculator error: " + str(e))
        */
    }

    double KmeansConcaveHull::haversineDistance(const lat_lon_coord first, const lat_lon_coord second)
    {
        const double earths_radius = 6371;

        // Get the difference between our two points then radians the difference into radians
        const double lat_delta = radians(second.Lat - first.Lat);
        const double lon_delta = radians(second.Lon - first.Lon);

        const double converted_lat1 = radians(first.Lat);
        const double converted_lat2 = radians(second.Lat);

        const double a =
            pow(sin(lat_delta / 2), 2) + cos(converted_lat1) * cos(converted_lat2) * pow(sin(lon_delta / 2), 2);

        const double c = 2 * atan2(sqrt(a), sqrt(1 - a));
        const double d = earths_radius * c;

        return d;
    }

    std::vector<double> KmeansConcaveHull::calculateDistances(lat_lon_coord currentPoint, const std::vector<lat_lon_coord>& kNearestPoints)
    {
        std::vector<double> distances;
        distances.reserve(kNearestPoints.size());
        for(int i = 0; i < kNearestPoints.size(); i++)
        {
           distances.push_back(haversineDistance(currentPoint, kNearestPoints[i]));
        }

        return distances;
    }

    uint32_t KmeansConcaveHull::getLowestLatitudeIndex()
    {
        std::vector<double> temp_lats = _lat;
        std::vector<double>::iterator it = std::min_element(std::begin(temp_lats), std::end(temp_lats));
        uint32_t index = std::distance(std::begin(temp_lats), it);
        return index;
    }

    std::vector<uint32_t> KmeansConcaveHull::getMaskedIndices(const std::vector<uint32_t>& input_array, const std::vector<bool>& mask)
    {
        std::vector<uint32_t> masked_array;
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
    std::vector<uint32_t> range(size_t size)
    {
        std::vector<uint32_t> output_array(size, 0);
        std::iota(output_array.begin(), output_array.end(), 0);
        return output_array;
    }

    std::vector<uint32_t> KmeansConcaveHull::getKNearest(uint32_t currentPointIndex, size_t k)
    {
        printf("CURRENT POINT INDEX: %u\n", currentPointIndex);
        printf("IX: [%.8f, %.8f]\n", _data_set[currentPointIndex].Lat, _data_set[currentPointIndex].Lon);

        _mask[6] = false;
        printf("IXS: [");
        for(int i = 0; i < _mask.size(); i++)
        {
            std::cout << _mask[i] << ", ";
        }
        std::cout << "]" << std::endl;

        std::vector<uint32_t> base_indices = getMaskedIndices(range(_mask.size()), _mask);
        printf("BASE_INDICES: [");
        for(auto& index: base_indices)
        {
            printf("%u, ", index);
        }
        std::cout << "]" << std::endl;

        printf("DATASET: [");
        for(int i = 0; i < _data_set.size(); i++)
        {
            printf("[%.8lf, %.8lf]\n", _data_set[i].Lat, _data_set[i].Lon);
        }

        std::vector<lat_lon_coord> masked_data_set = arraySubset(_data_set, base_indices);
        printf("MASKED_DATASET: [");
        for(auto& coord: masked_data_set)
        {
            printf("[%.8lf, %.8lf]\n", coord.Lat, coord.Lon);
        }
        //We need to calculate the distances array from the point to all other points
        //std::cout << masked_data_set[currentPointIndex].Lat << ", " << masked_data_set[currentPointIndex].Lon << std::endl;

        std::vector<double> distances = calculateDistances(masked_data_set[currentPointIndex], masked_data_set);

        for(auto& distance: distances)
        {
            std::cout << distance << std::endl;
        }

        //Sort the distances array in non-decending order
        std::vector<uint32_t> sorted_indices = argsort(distances);

        //Get the index of the lowest K points
        size_t k_check = std::min(_k, sorted_indices.size());
        sorted_indices.resize(k_check);

        //Set the index of the lowest K points to True, rest are false
        std::vector<uint32_t> kNearest = arraySubset(base_indices, sorted_indices);
        
        return kNearest;
    }

    template<typename T, typename K>
    std::vector<K> KmeansConcaveHull::arraySubset(const std::vector<T>& input_array, const std::vector<uint32_t>& indexes)
    {
        std::vector<K> output_array(indexes.size());

        for(int i = 0; i < indexes.size(); i++)
        {
            output_array[i] = input_array[indexes[i]];
        }

        return output_array;
    }

    /**
     * Argsort(currently support ascending sort)
     * @tparam T array element type
     * @param array input array
     * @return indices w.r.t sorted array
     */
    template<typename T>
    std::vector<uint32_t> KmeansConcaveHull::argsort(const std::vector<T> &array) {
        std::vector<uint32_t> indices(array.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(),
                [&array](int left, int right) -> bool {
                    // sort indices according to corresponding array element
                    return array[left] < array[right];
                });

        return indices;
    }

    double KmeansConcaveHull::getNextK()
    {
        if (_prime_ix < _prime_k.size())
        {
            return _prime_k[_prime_ix];
        }
        else
        {
            return -1.0f;
        }
    }

    std::vector<double> calculateHeadings(lat_lon_coord currentPointIndex, 
                                        const std::vector<lat_lon_coord>& searchPoints, 
                                        double ref_heading = 0.0l)
    {
        /*
        if ((ref_heading < 0.0f) || (ref_heading >= 360.0f))
        {
            //raise ValueError("The reference heading must be in the range [0, 360]")
        }
        
        float r_ix = radians(self.data_set[ix, :]);
        float r_ixs = radians(self.data_set[ixs, :]);

        float lon_dif = r_ixs[:, 0] - r_ix[0];

        float y = np.multiply(sin(lon_dif), cos(r_ixs[:, 1]));

        float x = cos(r_ix[1]) * sin(r_ixs[:, 1]) - sin(r_ix[1]) *
                np.multiply(cos(r_ixs[:, 1]), cos(lon_dif));

        float bearings = (degrees(arctan2(y, x)) + 360.0f) % 360.0f - ref_heading;

        bearings[bearings < 0.0] += 360.0;
        return bearings
        */
        return std::vector<double>(3, 1.0f);
    }

    bool containedCheck(const std::vector<lat_lon_coord>& hull, lat_lon_coord point)
    {
        return true;
    }
}