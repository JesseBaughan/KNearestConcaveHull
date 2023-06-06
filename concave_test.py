import os
import sys
import pandas as pd
import numpy as np
from descartes import PolygonPatch
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import alphashape
from TDSConcaveHull import TDSConcaveHull
import random as random
import cProfile, pstats

def plot_result(points_2d, alpha_shape):
    fig, ax = plt.subplots()
    ax.add_patch(alpha_shape)
    ax.scatter(*zip(*points_2d), color='red', zorder=1)
    plt.show()

def plot_result_alpha(points_2d, alpha_shape):
    fig, ax = plt.subplots()
    ax.add_patch(alpha_shape)
    ax.scatter(*zip(*points_2d), color='red', zorder=1)
    plt.show()

def plot_2d_points(points_2d):
    fig, ax = plt.subplots()
    ax.scatter(*zip(*points_2d))
    plt.show()

def generate_point(mean_x, mean_y, deviation_x, deviation_y):
    return random.gauss(mean_x, deviation_x), random.gauss(mean_y, deviation_y)

def generate_random_cluster():
    cluster_mean_x = 100
    cluster_mean_y = 100
    cluster_deviation_x = 50
    cluster_deviation_y = 50
    point_deviation_x = 5
    point_deviation_y = 5

    number_of_clusters = 1
    points_per_cluster = 100

    cluster_centers = [generate_point(cluster_mean_x,
                                    cluster_mean_y,
                                    cluster_deviation_x,
                                    cluster_deviation_y)
                    for i in range(number_of_clusters)]

    points = [generate_point(center_x,
                            center_y,
                            point_deviation_x,
                            point_deviation_y)
            for center_x, center_y in cluster_centers
            for i in range(points_per_cluster)]

    return points


def main():
    # Create a test cluster
    #points_2d = [[0., 0.], [0., 1.], [1., 1.], [1., 0.],
    #        [0.5, 0.25], [0.5, 0.75], [0.25, 0.5], [0.75, 0.5]]

    # Create a test cluster
    points_2d = [[-27.4507974, 153.0476219],
                [-27.451273651267112, 153.04687382688317], 
                [-27.45169316203486, 153.0476602627239], 
                [-27.451034639346872, 153.04883032580403],
                [-27.452475086130043, 153.05000678267146], 
                [-27.452644945552006, 153.04604902831844], 
                [-27.450159439942457, 153.04581245818747], 
                [-27.44996617171483, 153.04982136283908]]

    points_2d = np.array(points_2d)
    #points_2d = np.array(list(zip(*points_2d)))

    n_samples = 1
    failed_counter = 0
    tds_concave_hull = TDSConcaveHull() 

    profiler = cProfile.Profile()
    profiler.enable()

    for i in range(n_samples):
        #points_2d = np.array(generate_random_cluster())

        boundary = tds_concave_hull.calculate(points_2d)

        if boundary is None:
            #print("Failed")
            failed_counter += 1
        else:
            #shape_plot = plt.Polygon(boundary, fill = False, hatch='/')
            shape = TDSConcaveHull.boundary_points_to_shape(boundary)
            #plot_result(points_2d, shape_plot)

    print("FAILS: ", failed_counter)

    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats('cumtime')
    #stats.print_stats()

if __name__ == "__main__":
    main()