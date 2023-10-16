/**
 * @file ArrayUtils.cpp
 *
 * @brief Arrays utility functions used by concave hull calculation.
 * 
 * @author Jesse B
 */

#include "ArrayUtils.h"

using namespace std;

namespace Clustering
{

// Returns a range from 0 to size
vector<uint32_t> range(size_t size)
{
    vector<uint32_t> output_array(size, 0);
    iota(output_array.begin(), output_array.end(), 0);
    return output_array;
}

}

/*** end of file ***/