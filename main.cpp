#include <iostream>

#include "KmeansConcaveHull/KMeansConcaveHull.hpp"

int main()
{
    std::vector<float> lon = {-27.4507974, -27.451273651267112, -27.45169316203486, 
                              -27.451034639346872, -27.452475086130043, -27.452644945552006, 
                              -27.450159439942457, -27.44996617171483};

    std::vector<float> lat = {153.0476219, 153.04687382688317, 153.0476602627239, 
                              153.04883032580403, 153.05000678267146, 153.04604902831844, 
                              153.04581245818747, 153.04982136283908};

    Clustering::KmeansConcaveHull hull_calculator(lat, lon);

    const std::vector<bool> indices = hull_calculator.get_mask();
    for(int i = 0; i < indices.size(); i++)
    {
        std::cout << indices[i] << std::endl;
    }

    uint32_t lowest_lat_index = hull_calculator.getLowestLatitudeIndex();
    std::cout << "lowest lat index: " << lowest_lat_index << std::endl;


    return 0;
}