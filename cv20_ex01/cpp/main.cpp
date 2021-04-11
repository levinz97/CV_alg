#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <eigen3/Eigen/Core>
using namespace std;
using namespace cv;

int main() {
    double hx = 1;// grid size in both direction
    double hy = 1;
    double fx = 0;// gradient in both direction
    double fy = 0;
    double threshold = 0;// absolute value [0,255]

    std::cout << "OpenCV version:" <<CV_VERSION<< std::endl;
    string filename = "coins.pgm";
    Mat inimage,gradient_map;
    if(filename.substr(filename.find('.') + 1) == "pgm"){
        // read inimage with only grey scale
        inimage = cv::imread(filename, cv::IMREAD_GRAYSCALE);
        gradient_map = inimage.clone();
    }else{
        // for color image we need joint color gradient
        // TODO
        printf("unsupported type");
        return -1;
//        inimage = cv::imread(filename, cv::IMREAD_COLOR);
    }

    // create a dummy boundary by mirroring
    Mat new_image = Mat::zeros(inimage.rows+2, inimage.cols+2, inimage.type());
    inimage.copyTo(new_image(Rect(1,1,inimage.cols,inimage.rows)));
//    for (size_t x = 0; x < new_image.rows; x++){
//        new_image.at<uchar>(0,x) = new_image.at<uchar>(1,x);
//        new_image.at<uchar>(new_image.cols-1, x) = new_image.at<uchar>(new_image.cols-2, x);
//    }
//    for (size_t y = 0; y < new_image.cols; y++){
//        new_image.at<uchar>(y,0) = new_image.at<uchar>(y,1);
//        new_image.at<uchar>(y,new_image.rows-1) = new_image.at<uchar>(y, new_image.rows-2);
//    }
    inimage.row(0).copyTo(
            new_image(Rect(1,0,inimage.cols,1)));
    inimage.row(inimage.rows-1).copyTo(
            new_image(Rect(1,new_image.rows-1, inimage.cols,1)));
    inimage.col(0).copyTo(
            new_image(Rect(0,1,1, inimage.rows)));
    inimage.col(inimage.cols-1).copyTo(
            new_image(Rect(new_image.cols-1,1,1,inimage.rows)));
    // trivial, mirror 4 corners, needed for 2D filters
    new_image.at<uchar>(0,0) = inimage.at<uchar>(0,0);
    new_image.at<uchar>(0,new_image.cols-1) = inimage.at<uchar>(0,inimage.cols-1);
    new_image.at<uchar>(new_image.rows-1,0) = inimage.at<uchar>(inimage.rows-1,0);
    new_image.at<uchar>(new_image.rows-1,new_image.cols-1) = inimage.at<uchar>(inimage.rows-1,inimage.cols-1);

    Mat sm_image = new_image.clone();
    // Gaussian pre-smoothing note that size value must be odd
    GaussianBlur(new_image,sm_image,Size(3,3),200);
    imshow("original image", inimage);
    imshow("pre-smoothed image", sm_image);
    waitKey();

    //compute the gradiant_filter
    Eigen::Vector3d gradiant_filter = {-1, 0, 1};
    gradiant_filter /= 2;
    Eigen::Vector3d Intensity = Eigen::Vector3d::Zero();
//    Intensity << new_image.at<uchar>(1,0), new_image.at<uchar>(1,1),new_image.at<uchar>(1,2);
    for(size_t y = 1; y < sm_image.rows - 2; y++){
        for(size_t x = 1; x < sm_image.cols - 2; x++){
            // compute gradient in x direction
            for(int t = -1; t < Intensity.size() - 1; t++){
                Intensity[t + 1] = sm_image.at<uchar>(y, x+t);
            }
            fx = Intensity.dot(gradiant_filter / hx);
            // compute gradient in y direction
            for(int t = -1; t < Intensity.size() - 1; t++){
                Intensity[t + 1] = sm_image.at<uchar>(y + t, x);
            }
            fy = Intensity.dot(gradiant_filter / hy);
            // thresholding the magnitude of gradient, set the minors to 0
            double magnitude =sqrt(pow(fx, 2) + pow(fy, 2));
            if(magnitude < threshold) magnitude = 0;
            gradient_map.at<uchar>(y - 1,x - 1) = magnitude;

        }
    }

//    gradient_map(Rect(0,0,100,100)).setTo(Scalar(125));
//    cout<<(int) gradient_map.at<uchar>(-1,1) <<endl;
    imshow("gradient_map", gradient_map);
//    imshow("new_image",new_image);
    waitKey(0);

    return 0;
}
