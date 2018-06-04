#include "Slic.h"

/*Constructor. Nothing is done here*/
Slic::Slic()
{

}

/*Destructor. Clear any present data*/
Slic::~Slic()
{
	clear_data();
}

/*Clear the data as saved by the algorithm

* Input: - 
* Output: -
*/
void Slic::clear_data()
{
	clusters.clear();
	distances.clear();
	centers.clear();
	center_counts.clear();
}

/*Initialize the cluster centers and initial values of the pixle-wise cluster*/
void Slic::init_data(Mat image)
{
	int hei = image.rows;
	int wid = image.cols;
	/*Initialize the cluster and distance matrices.*/
	for (int i = 0; i < wid; i++)
	{
		vector<int> cr;
		vector<double> dr;
		for (int j = 0; j < hei; j++)
		{
			cr.push_back(-1);
			dr.push_back(FLT_MAX);
		}
		clusters.push_back(cr);
		distances.push_back(dr);
	}
	/*Initialize the centers and counters.*/
	for (int i = step / 2; i < wid - step / 2; i += step)
	{
		for(int j = step / 2; j < hei - step / 2; j += step)
		{
			vector<double> center;
			/*Find the local minimum (gradient-wise)*/
			Point nc = find_local_minimum(image, Point(i, j));
			Vec3b color = image.at<Vec3b>(nc.y, nc.x);

			/*Generate the center vector*/
			center.push_back(color[0]);
			center.push_back(color[1]);
			center.push_back(color[2]);
			center.push_back(nc.x);
			center.push_back(nc.y);

			/*Append to vector of centers.*/
			centers.push_back(center);
			center_counts.push_back(0);
		}
	}

}

/*Compute the distance between a cluser center and an individual pixel.

Input: The cluster index(int), the pixel(Point), and the Lab values of the pixel
Output: The distances(double)
*/
double Slic::compute_dist(int ci, Point pixel, Vec3b color)
{
	double dc = sqrt(pow(centers[ci][0] - color[0], 2) + pow(centers[ci][1] - color[1], 2)
		             + pow(centers[ci][2] - color[2], 2));
	double ds = sqrt(pow(centers[ci][3] - pixel.x, 2) + pow(centers[ci][4] - pixel.y, 2));

	return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
}


/*Fine a local gradient minimum of a pixel in a 3x3 neighborhood .
This method is called upon initialization of the cluster centers.

Input: The image and the pixel center
Output: The local gradient minimum.
*/
Point Slic::find_local_minimum(Mat image, Point center)
{
	double min_grad = FLT_MAX;
	Point loc_min = Point(center.x, center.y);

	for (int i = center.x - 1; i < center.x + 2; i++)//point.x = cols
	{

		for (int j = center.y - 1; j < center.y + 2; j++)
		{
			Vec3b c1 = image.at<Vec3b>(j+1, i);
			Vec3b c2 = image.at<Vec3b>(j, i + 1);
			Vec3b c3 = image.at<Vec3b>(j, i);
			/*Convert color value to grayscale values.*/
			double i1 = c1[0] * 0.114 + c1[1] * 0.587 + c1[2] * 0.299;
			double i2 = c2[0] * 0.114 + c2[1] * 0.587 + c2[2] * 0.299;
			double i3 = c3[0] * 0.114 + c3[1] * 0.587 + c3[2] * 0.299;

			/*Compute horizontal and vertical gradients and keep track 
			of the minimum*/
			if (sqrt(pow(i1 - i3, 2) + pow(i2 - i3, 2)) < min_grad)
			{
				//min_grad = fabs(i1 - i3) + fabs(i2 - i3);
				min_grad = sqrt(pow(i1 - i3, 2) + pow(i2 - i3, 2));
				loc_min.x = i;
				loc_min.y = j;
			}
		}
	}
	return loc_min;
}

/*Compute the over-segmentation based on the step-size and relative weighting
 of the pixel and color values
 
 Input: The Lab image, the stepSize, and the weight(int)
 Output: -
 */
void Slic::generate_superpixels(Mat image, int step, int nc)
{
	this->step = step;
	this->nc = nc;
	this->ns = step;

	int hei = image.rows;
	int wid = image.cols;

	/*Clear previous data(if any), and re-initialize it.*/
	clear_data();
	init_data(image);

	/*Run EM for 10 iterations (as prescribed by the algorithm).*/
	for (int k = 0; k < NR_ITERATIONS; k++)
	{
		/*Reset distance values.*/
		for (int j = 0; j < wid; j++)
		{
			for (int i = 0; i < hei; i++)
			{
				distances[j][i] = FLT_MAX;
			}
		}

		for (int ci = 0; ci < (int)centers.size(); ci++)
		{
			/*Only compare to pixels in a 2 x step by 2 x step region.*/
			for (int j = centers[ci][3] - step; j < centers[ci][3] + step + 1; j++)
			{
				for (int i = centers[ci][4] - step; i < centers[ci][4] + step + 1; i++)
				{
					if (j >= 0 && j < wid && i >= 0 && i < hei)
					{
						Vec3b color = image.at<Vec3b>(i, j);
						double d = compute_dist(ci, Point(j, i), color);

						/*Update cluster allocation if the cluster minimizes the distance*/
						if (d < distances[j][i])//可以将距离换成一个恒定的值
						{
							distances[j][i] = d;
							clusters[j][i] = ci;
						}
					}
				}
			}
		}

		/*Clear the center values*/
		for (int ci = 0; ci < (int)centers.size(); ci++)
		{
			for (int i = 0; i < 5; i++)
			{
				centers[ci][i] = 0;
			}
			center_counts[ci] = 0;
		}

		/*Compute the new cluster centers.*/
		for (int j = 0; j < wid; j++)
		{
			for (int i = 0; i < hei; i++)
			{
				int c_id = clusters[j][i];

				if (c_id != -1)
				{
					Vec3b color = image.at<Vec3b>(i, j);

					centers[c_id][0] += color[0];
					centers[c_id][1] += color[1];
					centers[c_id][2] += color[2];
					centers[c_id][3] += j;
					centers[c_id][4] += i;

					center_counts[c_id] += 1;
				}
			}
		}
		/*Normalize the clusters.*/
		for (int ci = 0; ci < (int)centers.size(); ci++)
		{
			centers[ci][0] /= center_counts[ci];
			centers[ci][1] /= center_counts[ci];
			centers[ci][2] /= center_counts[ci];
			centers[ci][3] /= center_counts[ci];
			centers[ci][4] /= center_counts[ci];
		}
	}
}

/*Enforce connectivity of the superpixels.
This part is not actively discussed in the paper, but forms an active part of the implementation of the 
authors of the paper

Input: The image
Output: -
*/
void Slic::create_connectivity(Mat image)
{
	int label = 0, adjlabel = 0;
	int hei = image.rows;
	int wid = image.cols;
	const int lims = (hei * wid) / ((int)centers.size());

	const int dx4[4] = { -1, 0, 1, 0 };
	const int dy4[4] = { 0, -1, 0, 1 };

	/*Initialize the new cluster matrix*/
	vec2di new_clusters;
	for (int j = 0; j < wid; j++)
	{
		vector<int> nc;
		for (int i = 0; i < hei; i++)
		{
			nc.push_back(-1);
		}
		new_clusters.push_back(nc);
	}
	for (int j = 0; j < wid; j++)
	{
		for (int i = 0; i < hei; i++)
		{
			if (new_clusters[j][i] = -1)
			{
				vector<Point> elements;
				elements.push_back(Point(j, i));

				/*Find an adjacent label. for possible use later.*/
				for (int k = 0; k < 4; k++)
				{
					int x = elements[0].x + dx4[k], y = elements[0].y + dy4[k];

					if (x >= 0 && x < wid && y >= 0 && y < hei)
					{
						if (new_clusters[x][y] >= 0)
						{
							adjlabel = new_clusters[x][y];
						}
					}
				}

				int count = 1;
				for (int c = 0; c < count; c++)
				{
					for (int k = 0; k < 4; k++)
					{
						int x = elements[c].x + dx4[k], y = elements[c].y + dy4[k];
						if (x >= 0 && x < wid && y >= 0 && y < hei)
						{
							if (new_clusters[x][y] == -1 && clusters[j][i] == clusters[x][y])
							{
								elements.push_back(Point(x, y));
								new_clusters[x][y] = label;
								count += 1;
							}
						}
					}
				}

				/*Use the earlier found adjacent label if a segment size is smaller than a limit*/
				if (count <= lims >> 2)
				{
					for (int c = 0; c < count; c++)
					{
						new_clusters[elements[c].x][elements[c].y] = adjlabel;
					}
					label -= 1;
				}
				label += 1;
			}
		}
	}
}

/*Display a single pixel wide contour around the clusters.
Input: The target image and contour color
Output: -
*/
void Slic::display_contours(Mat image, Scalar color)
{
	const int dx8[8]{-1, -1, 0, 1, 1, 1, 0, -1 };
	const int dy8[8]{0, -1, -1, -1, 0, 1, 1, 1 };
	int hei = image.rows;
	int wid = image.cols;
	
	/*Initialize the contour vector and the matrix detailing whether 
	a pixel is already taken to be a contour*/
	vector<Point> contours;
	vec2db istaken;
	for (int j = 0; j < wid; j++)
	{
		vector<bool> nb;
		for (int i = 0; i < hei; i++)
		{
			nb.push_back(false);
		}
		istaken.push_back(nb);
	}

	/*Go through all the pixels*/
	for (int j = 0; j < wid; j++)
	{
		for (int i = 0; i < hei; i++)
		{
			int nr_p = 0;

			/*Compare the pixel to its 8 neighbors.*/
			for (int k = 0; k < 8; k++)
			{
				int x = j + dx8[k], y = i + dy8[k];

				if (x >= 0 && x < wid && y >= 0 && y < hei)
				{
					if (istaken[x][y] == false && clusters[j][i] != clusters[x][y])
					{
						nr_p += 1;
					}
				}
			}
			/*Add the pixel to the contour list if desired*/
			if (nr_p >= 3)
			{
				contours.push_back(Point(j, i));
				istaken[j][i] = true;
			}
		}
	}
	/*Draw the contour pixels*/
	for (int i = 0; i < (int)contours.size(); i++)
	{
		int x = contours[i].x;
		int y = contours[i].y;
		image.at<Vec3b>(y, x) = Vec3b(255, 0, 0);
		//circle(image, contours[i], 1, color);
		//imshow("image", image);
		//waitKey(100);

	}
}

/*Given the pixels of each cluster the same color values.
The specified color is the mean RGB color per cluster.

Input: The target iamge
Output: -
*/
void Slic::color_with_cluster_means(Mat image)
{
	vector<Scalar> colors(centers.size());
	int hei = image.rows;
	int wid = image.cols;

	/*Gather the color values per cluster.*/
	for (int j = 0; j < wid; j++)
	{
		for (int i = 0; i < hei; i++)
		{
			int index = clusters[j][i];
			Vec3b color = image.at<Vec3b>(i, j);

			colors[index][0] += color[0];
			colors[index][1] += color[1];
			colors[index][2] += color[2];
		}
	}

	/*Divide by the number of pixels per cluster to get the mean color*/
	for (int i = 0; i < (int)colors.size(); i++)
	{
		colors[i][0] /= center_counts[i];
		colors[i][1] /= center_counts[i];
		colors[i][2] /= center_counts[i];
	}
	
	/*Fill in.*/
	for (int j = 0; j < wid; j++)
	{
		for (int i = 0; i < hei; i++)
		{
			Scalar ncolor = colors[clusters[j][i]];
			circle(image, Point(j, i), 1, ncolor);
		}
	}
}