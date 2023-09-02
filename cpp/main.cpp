#include <iostream>

#include "KMeansConcaveHull/KMeansConcaveHull.hpp" 


int main()
{
    std::vector<double> lat = {-27.4507974, -27.451273651267112, -27.45169316203486, 
                                -27.451034639346872, -27.452475086130043, -27.452644945552006, 
                                -27.450159439942457, -27.44996617171483};

    std::vector<double> lon = {153.0476219, 153.04687382688317, 153.0476602627239, 
                                153.04883032580403, 153.05000678267146, 153.04604902831844, 
                                153.04581245818747, 153.04982136283908};

    Clustering::KmeansConcaveHull hullCalc(lat, lon);

    uint32_t index = hullCalc.getLowestLatitudeIndex();

    std::vector<uint32_t> knn = hullCalc.getKNearest(index);

    std::vector<double> headings = hullCalc.calculateHeadings(index, knn, 270.0l);
    Clustering::NegateArray<double>(headings);

    for(const auto& heading : headings)
    {
        std::cout << "Value: " << heading << std::endl;
    }

    hullCalc.calculate(3);

	Clustering::pdd A(1, 1);
	Clustering::pdd B(4, 4);
	Clustering::pdd C(1, 7);
	Clustering::pdd D(3, 4);
    Clustering::pdd intersectPoint = Clustering::lineLineIntersection(A, B, C, D);

	if (intersectPoint.first == FLT_MAX &&
		intersectPoint.second==FLT_MAX)
	{
		std::cout << "The given lines AB and CD are parallel.\n";
	}

	else
	{
		// NOTE: Further check can be applied in case
		// of line segments. Here, we have considered AB
		// and CD as lines
		std::cout << "The intersection of the given lines AB "
				"and CD is: ";
	}

    return 0;
}