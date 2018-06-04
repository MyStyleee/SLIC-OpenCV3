//#pragma once
#ifndef SLIC_H
#define SLIC_H

#include "CommonFunc.h"

/*2d matrices are handled by 2d vectors*/
#define vec2dd vector<vector<double>>
#define vec2di vector<vector<int>>
#define vec2db vector<vector<bool>>
/*The number of iterations run by the clustring algorithm*/
#define NR_ITERATIONS 10

class Slic {
	private:
		/*The cluster assignments and distance values for each pixel*/
		vec2di clusters;
		vec2dd distances;

		/*The LAB and xy values of the centers.*/
		vec2dd centers;
		/*The number of occurences of each center.*/
		vector<int> center_counts;

		/*The step size per cluster, and the color(nc) and distance(ns) parameters*/
		int step, nc, ns;

		/*Remove and Initialize the 2d vectors.*/
		void clear_data();
		void init_data(Mat image);

		/*Compute the distance between a center and an individual pixel*/
		double compute_dist(int ci, Point pixel, Vec3b color);
		/*Find the pixel with the lowest gradient in a 3x3 surroundings.*/
		Point find_local_minimum(Mat image, Point center);

	public:
		/*Class constructors and deconstructors*/
		Slic();
		~Slic();

		/*Generate an over-segmentation for an image.*/
		void generate_superpixels(Mat image, int step, int nc);
		/*Enforce connectivity for an image.*/
		void create_connectivity(Mat image);

		/*Draw functions, Resp, dispalyal of the centers and the contours.*/
		void display_center_grid(Mat image, Scalar color);
		void display_contours(Mat image, Scalar color);
		void color_with_cluster_means(Mat image);

};
#endif