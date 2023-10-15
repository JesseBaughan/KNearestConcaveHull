#include "KMeansConcaveHull/KMeansConcaveHull.hpp" 

int main()
{
    vector<double> lat = {-27.4507974, -27.451273651267112, -27.45169316203486, 
                                -27.451034639346872, -27.452475086130043, -27.452644945552006, 
                                -27.450159439942457, -27.44996617171483};

    vector<double> lon = {153.0476219, 153.04687382688317, 153.0476602627239, 
                                153.04883032580403, 153.05000678267146, 153.04604902831844, 
                                153.04581245818747, 153.04982136283908};

    // Example usage
    Clustering::KmeansConcaveHull hullCalc(lat, lon);
    vector<Clustering::Point> hull = hullCalc.calculate();

    return 0;
}