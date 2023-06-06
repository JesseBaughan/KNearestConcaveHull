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
#include <algorithm>
#include <numeric>

#include "KMeansConcaveHull.hpp"

namespace Clustering
{
    static inline float radians(float degrees)
    {
        return (degrees * M_PI)  / 180.0f;
    }

    static inline float degrees(float radians)
    {
        return radians * (180.0 / M_PI);
    }

    std::vector<std::vector<float>> KmeansConcaveHull::calculate(const std::vector<std::vector<float>> &points, uint32_t k)
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

    float KmeansConcaveHull::haversineDistance(const lat_lon_coord first, const lat_lon_coord second)
    {
        const float earths_radius = 6371;

        // Get the difference between our two points then radians the difference into radians
        const float lat_delta = radians(second.Lat - first.Lat);
        const float lon_delta = radians(second.Lon - first.Lon);

        const float converted_lat1 = radians(first.Lat);
        const float converted_lat2 = radians(second.Lat);

        const float a =
            pow(sin(lat_delta / 2), 2) + cos(converted_lat1) * cos(converted_lat2) * pow(sin(lon_delta / 2), 2);

        const float c = 2 * atan2(sqrt(a), sqrt(1 - a));
        const float d = earths_radius * c;

        return d;
    }

    std::vector<float> KmeansConcaveHull::calculateDistances(lat_lon_coord currentPoint, const std::vector<lat_lon_coord>& kNearestPoints)
    {
        std::vector<float> distances;
        distances.reserve(kNearestPoints.size());
        for(int i = 0; i < kNearestPoints.size(); i++)
        {
           distances.push_back(haversineDistance(currentPoint, kNearestPoints[i]));
        }

        return distances;
    }

    uint32_t KmeansConcaveHull::getLowestLatitudeIndex()
    {
        std::vector<float> temp_lats = _lat;
        std::vector<float>::iterator it = std::min_element(std::begin(temp_lats), std::end(temp_lats));
        uint32_t index = std::distance(std::begin(temp_lats), it);
        return index;
    }

    template<typename Type>
    std::vector<Type> KmeansConcaveHull::getMaskedArray(const std::vector<Type>& input_array, std::vector<bool>& mask)
    {
        std::vector<Type> masked_array;
        masked_array.reserve(input_array.size());

        for(int i = 0; i < array.size(); i++)
        {
            if(mask[i] == true)
            {
                masked_array.push_back(increasing_vals[i]);
            }
        }

        return masked_array;
    }

    std::vector<bool> KmeansConcaveHull::getKNearest(uint32_t currentPointIndex, uint32_t k)
    {
        std::vector<bool> ixs = _indices;
        std::vector<uint32_t> increasing_vals(ixs.size(), 0);
        std::iota(increasing_vals.begin(), increasing_vals.end(), 0);

        std::vector<uint32_t> base_indices = getMaskedArray(increasing_vals, ixs);

        std::vector<lat_lon_coord> masked_data_set;
        masked_data_set.reserve(base_indices.size());
        for(int i = 0; i < base_indices.size(); i++)
        {
            masked_data_set[i] = _data_set[base_indices[i]];
        }

        std::vector<float> distances = calculateDistances(masked_data_set[currentPointIndex], masked_data_set);

        sorted_indices = np.argsort(distances)

        k_check = min(k, len(sorted_indices))
        k_nearest = sorted_indices[range(k_check)]
        return base_indices[k_nearest]
        return std::vector<bool>(3, true);
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

    std::vector<float> calculateHeadings(lat_lon_coord currentPointIndex, 
                                        const std::vector<lat_lon_coord>& searchPoints, 
                                        float ref_heading=0.0f)
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
    return std::vector<float>(3, 1.0f);
    }

    bool containedCheck(const std::vector<lat_lon_coord>& hull, lat_lon_coord point)
    {
        return true;
    }
}