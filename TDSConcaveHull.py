import math
import numpy as np

from shapely.geometry import Point, LineString, Polygon


"""
@class TDSConcaveHull

@brief A class to calculate the Concave hull of a set of points.
"""
class TDSConcaveHull():

    """
    @fn __init__

    @brief The initialisation function for the TDSConcaveHull class.

    @param points The set of points to find the concave hull of.
    @param prime_ix The index of the prime array to start with.
    """
    def __init__(self, prime_ix=0):
        
        # prime k's which have been observed being used
        # 43, 53, 61, 71, 79, 89, 

        self.prime_k = np.array([3,  7,   13,  19,  29,  37,
                                 43, 53,  61,  71,  79,  89,
                                 97, 101, 107, 113, 131, 139])
        
        self.prime_ix = prime_ix


    """
    @fn get_next_k

    @brief Returns the next k value in the prime series
    """
    def get_next_k(self):
        if self.prime_ix < len(self.prime_k):
            return self.prime_k[self.prime_ix]
        else:
            return -1

    """
    @fn haversine_distance

    @brief Calculates the haversine distance between two geological points.

    @param loc_ini The first geological point of shape [lon, lat].
    @param loc_end The sencond geological point of shape [lon, lat].
    """
    def haversine_distance(self, loc_ini, loc_end):
        lon1, lat1, lon2, lat2, = map(np.radians, [loc_ini[0], loc_ini[1], 
                                                   loc_end[:, 0], loc_end[:, 1]])
        
        lon_dif = abs(lon2 - lon1)
        lat_dif = abs(lat2 - lat1)

        a = np.square(np.sin(lat_dif / 2.0)) + np.cos(lat1) * np.cos(lat2) * \
            np.square(np.sin(lon_dif / 2.0))
        
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

        return 6371000.0 * c


    """
    @fn get_lowest_latitude_index

    @brief Gets the point with the lowest latitude.

    @param points The array of points to find the lowest point of.
    """
    @staticmethod
    def get_lowest_latitude_index(points):
        return np.argsort(points[:, 1])[0]

    """
    @fn get_k_nearest

    @brief Gets the k nearest neighbours

    @param ix The index of the node to search around.
    @param k The amount of neighbours to search around
    """
    def get_k_nearest(self, ix, k):
        print("CURRENT POINT INDEX", ix)
        print("K: ", k)

        ixs = self.indices

        arang = np.arange(len(ixs))
        base_indices = arang[ixs]

        print("IX: ", self.data_set[ix, :])
        print("MASKED DATASET: ", self.data_set[ixs, :])
        distances = self.haversine_distance(self.data_set[ix, :], self.data_set[ixs, :])
        sorted_indices = np.argsort(distances)
        print("DISTANCES: ", distances)
        print("ARGSORTED: ", sorted_indices)

        k_check = min(k, len(sorted_indices))
        k_nearest = sorted_indices[range(k_check)]
        print("KNEAREST: ", k_nearest)
        return base_indices[k_nearest]
    
    """
    @fn calculate_headings

    @brief Calculates the headings from a source point to a set of target points.

    @param ix Index to the source point in the dataset.
    @param ixs Indexes to calculate the bearing to from ix.
    @param ref_heading Reference heading measured in degrees counterclockwise from North.
    """
    def calculate_headings(self, ix, ixs, ref_heading=0.0):
        if ref_heading < 0 or ref_heading >= 360.0:
            raise ValueError("The reference heading must be in the range [0, 360]")
        
        r_ix = np.radians(self.data_set[ix, :])
        r_ixs = np.radians(self.data_set[ixs, :])

        lon_dif = r_ixs[:, 0] - r_ix[0]

        y = np.multiply(np.sin(lon_dif), np.cos(r_ixs[:, 1]))
        x = math.cos(r_ix[1]) * np.sin(r_ixs[:, 1]) - math.sin(r_ix[1]) * \
            np.multiply(np.cos(r_ixs[:, 1]), np.cos(lon_dif))
        bearings = (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0 - ref_heading

        bearings[bearings< 0.0] += 360.0
        return bearings

    """
    @fn recurse_calculate

    @brief Reruns the calculate another time with a larger amount of neighbours
    """
    def recurse_calculate(self):
        recurse = TDSConcaveHull(self.prime_ix + 1)
        next_k = recurse.get_next_k()
        if next_k == -1:
            #timing_logger.info("Hull calculator failed with size " + str(self.data_set.size))
            return None
        return recurse.calculate(self.data_set, next_k)

    """
    @fn calculate

    @brief Calculates the concave hull of the cluster.

    @param k The amount of neighbours to check against
    """
    def calculate(self, points, k=3):

        if isinstance(points, np.core.ndarray):
            self.data_set = points
        else:
            raise ValueError("Please provide a [N,2] numpy array")
        self.data_set = np.unique(self.data_set, axis=0)
       
        self.indices = np.ones(self.data_set.shape[0], dtype=bool)

        print(self.indices)

        try:
            if self.data_set.shape[0] <= 3:
                print("Skipped hull calc <= 3 points")
                #timing_logger.info("Skipped hull calc <= 3 points")
                return self.data_set
            
            k_check = min(k, self.data_set.shape[0]) 

            first_point = self.get_lowest_latitude_index(self.data_set)
            print("FIRST POINT: ", first_point)
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
                print("KNEAREST: ", knn)

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
                #timing_logger.info("Hull size: " + str(hull.size) + ", Prime K: " + str(self.prime_k[self.prime_ix]))
                return hull
            else: 
                return self.recurse_calculate()
        except Exception as e:
            print("HullCalculator error: " + str(e))


    def __contained_check(self, hull, point):
        poly = Polygon(hull)
        pt = Point(point)
        
        if poly.intersects(pt) or pt.within(poly):
            return True
        else:
            d = poly.distance(pt)
            if d < 1e-5:
                return True
            else:
                return False

    @staticmethod
    def boundary_points_to_shape(hull):
        if len(hull) > 2:
            return Polygon(hull)
        elif len(hull) == 2:
            return LineString(hull)
        else:
            return Point(hull)