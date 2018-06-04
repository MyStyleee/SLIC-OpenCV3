#include "Slic.h"

int main(int argc, char* argv[])
{
	double time = static_cast<double>(getTickCount());
	Mat image = imread(argv[1], 1);
	Mat lab_image = image.clone();
	cvtColor(image, lab_image, CV_BGR2Lab);

	cout << "begin****" << endl;
	/*Yield the number of superpixels and weight-factors from the user*/
	int h = image.rows, w = image.cols;
	int nr_superpixels = atoi(argv[2]);
	int nc = atoi(argv[3]);

	double step = sqrt((w*h) / (double)nr_superpixels);

	/*Perform the SLIC superpixel algorithm*/
	Slic slic;
	slic.generate_superpixels(lab_image, step, nc);
	slic.create_connectivity(lab_image);

	/*Fill the image with cluster means color*/
	slic.color_with_cluster_means(image);
	/*Display the contours and show the result*/
	slic.display_contours(image, Scalar(255, 0, 0));
	

	//circle(image, Point(100, 100), 1, Scalar(255, 0, 0));
	imshow("segment", image);
	waitKey(5);
	imwrite(argv[4], image);

	time = static_cast<double>(getTickCount()) - time;
	time /= getTickFrequency();
	cout << "total time is: " << time << endl;
	system("pause");
	return 0;
}