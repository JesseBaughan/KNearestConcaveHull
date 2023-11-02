#include "gtest/gtest.h"
#include "../KMeansConcaveHull/KMeansConcaveHull.h"


TEST(testCubeCluster, test1) 
{
    vector<double> lat = {-27.4507974, -27.451273651267112, -27.45169316203486, 
                                -27.451034639346872, -27.452475086130043, -27.452644945552006, 
                                -27.450159439942457, -27.44996617171483};

    vector<double> lon = {153.0476219, 153.04687382688317, 153.0476602627239, 
                                153.04883032580403, 153.05000678267146, 153.04604902831844, 
                                153.04581245818747, 153.04982136283908};

    Clustering::KmeansConcaveHull hullCalc(lat, lon);
    vector<Clustering::Point> hull = hullCalc.calculate(3);

    vector<Clustering::Point> answer = {{153.046049028318436, -27.4526449455520059}, 
                                        {153.05000678267146, -27.4524750861300433}, 
                                        {153.049821362839083, -27.4499661717148307}, 
                                        {153.04581245818747, -27.4501594399424569}, 
                                        {153.046873826883171, -27.4512736512671118}, 
                                        {153.047660262723895, -27.4516931620348608}, 
                                        {153.048830325804033, -27.4510346393468723}, 
                                        {153.047621899999996, -27.450797399999999}, 
                                        {153.046049028318436, -27.4526449455520059}};

    EXPECT_TRUE(hull.size() == answer.size());

    for(int i = 0; i < answer.size(); i++)
    {
        ASSERT_DOUBLE_EQ(hull[i].x, answer[i].x);
        ASSERT_DOUBLE_EQ(hull[i].y, answer[i].y);
    }
}

TEST(testStarCluster, test2) 
{
    vector<double> lat = {-29.3693, -29.4052, -29.4046, -29.43989, -29.48712, 
                          -29.45364, -29.49429, -29.44109, -29.41956, -29.407};

    vector<double> lon = {149.83226, 149.86041, 149.92839, 149.87963, 149.90435, 
                          149.83912, 149.78831, 149.79518, 149.73818, 149.81097};

    Clustering::KmeansConcaveHull hullCalc{lat, lon};
    vector<Clustering::Point> hull = hullCalc.calculate();

    vector<Clustering::Point> answer = {{149.78831, -29.49429},
                                        {149.90435, -29.48712},
                                        {149.92839, -29.4046 },
                                        {149.83226, -29.3693 },
                                        {149.73818, -29.41956},
                                        {149.79518, -29.44109},
                                        {149.83912, -29.45364},
                                        {149.87963, -29.43989},
                                        {149.78831, -29.49429},
                                        {149.78831, -29.49429}};

    EXPECT_TRUE(hull.size() == answer.size());

    for(int i = 0; i < answer.size(); i++)
    {
        ASSERT_DOUBLE_EQ(hull[i].x, answer[i].x);
        ASSERT_DOUBLE_EQ(hull[i].y, answer[i].y);
    }
}

 