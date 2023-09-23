#include <iostream>

#include "KMeansConcaveHull/KMeansConcaveHull.hpp" 
#include "LineIntersect.hpp"

int main()
{
    vector<double> lat = {-27.4507974, -27.451273651267112, -27.45169316203486, 
                                -27.451034639346872, -27.452475086130043, -27.452644945552006, 
                                -27.450159439942457, -27.44996617171483};

    vector<double> lon = {153.0476219, 153.04687382688317, 153.0476602627239, 
                                153.04883032580403, 153.05000678267146, 153.04604902831844, 
                                153.04581245818747, 153.04982136283908};

    Clustering::KmeansConcaveHull hullCalc(lat, lon);

    uint32_t index = hullCalc.getLowestLatitudeIndex(hullCalc._data_set);

    vector<uint32_t> knn = hullCalc.getKNearest(index);

    vector<double> headings = hullCalc.calculateHeadings(index, knn, 270.0l);
    Clustering::NegateArray<double>(headings);

    for(const auto& heading : headings)
    {
        cout << "Value: " << heading << endl;
    }

    hullCalc.calculate(hullCalc._data_set, 3);

	Point A(0, 1);
	Point B(2, 1);
	Point C(0, 3);
	Point D(1, 0);

    bool intersects = lineLineIntersection(A, B, C, D);
    if(intersects)
    {
        cout << "Intersects";
    }
    else
    {

        cout << "Does not intersect";
    }

    return 0;
}