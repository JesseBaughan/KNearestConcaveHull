/**
 * @file ArrayUtils.h
 *
 * @brief Arrays utility functions used by concave hull calculation.
 * 
 * @author Jesse B 
 */

#ifndef ARRAY_UTILS_H_
#define ARRAY_UTILS_H_

#include <stdint.h>
#include <vector>
#include <array>
#include <numeric>

namespace Clustering
{

std::vector<uint32_t> range(size_t size);

/**
 * Argsort(currently support ascending sort)
 * @tparam T array element type
 * @param array input array
 * @return indices w.r.t sorted array
 */
template<typename T>
std::vector<uint32_t> argsort(const std::vector<T> &array) {
    std::vector<uint32_t> indices(array.size());
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(),
            [&array](int left, int right) -> bool {
                return array[left] < array[right];
            });

    return indices;
}

template<typename T>
void NegateArray(std::vector<T>& inputVector)
{
    for(int i = 0; i < inputVector.size(); i++)
    {
        inputVector[i] = -inputVector[i];
    }
}

template<typename Type>
std::vector<Type> arraySubset(const std::vector<Type>& input_array, const std::vector<uint32_t>& indexes)
{
    std::vector<Type> output_array;
    output_array.reserve(indexes.size());
    for(int i = 0; i < indexes.size(); i++)
    {
        output_array.push_back(input_array[indexes[i]]);
    }

    return output_array;
}

}

#endif /* ARRAY_UTILS_H_ */

/*** end of file ***/