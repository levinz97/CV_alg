#include <iostream>
#include <chrono>
using namespace std;

#include <opencv2/core/core.hpp>

#include <opencv2/highgui/highgui.hpp>

int main() {
    float offset = 150;
    float scaling = 1.2;

    // read image with only grey scale
    cv::Mat image = cv::imread("head.pgm", cv::IMREAD_GRAYSCALE);

    if(image.data == nullptr){
        cerr<<"file doesn't exist"<<endl;
        return 0;
    }

    cout<<"breadth of image is "<< image.cols<<",image height is "<< image.rows <<" number of channel is "
    <<image.channels()<<endl;

    cv::imshow("simple_filter", image);
    cv::waitKey(0);//wait for keyboard input

    //clone the image for output, not directly " = "
    cv::Mat image_out = image.clone();

    //use std::chrono to count the time
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

    for (size_t y = 0; y < image_out.rows; y++){
        for (size_t x = 0; x < image_out.cols; x++){
            auto intensity =  (float)image_out.at<uchar>(y, x);
//            cout << intensity << endl;
            if (intensity + offset > 255){
                image_out.at<uchar>(y,x)= 255;
            }
            else{
                intensity + offset < 0 ? image_out.at<uchar>(y,x) = 0 : image_out.at<uchar>(y,x) += offset;
            }
        }
    }

    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>> (t2 - t1);
    cout<<"used time is "<< time_used.count() <<endl;
    cv::imshow("output", image_out);
    cv::waitKey();
    cv::imwrite("head_out.pgm", image_out);

    cv::Mat image_out1 = image.clone();
    for (size_t y = 0; y < image_out1.rows; y++){
        auto *row_ptr = image.ptr<float>(y);
        for(size_t x = 0; x < image_out1.cols; x++) {
            auto *data_ptr = &row_ptr[x];
            *data_ptr - 127.5 < 0 ? *data_ptr = 0 : *data_ptr -= 127.5;
            *data_ptr *= scaling * *data_ptr;
            *data_ptr + 127.5 > 255? *data_ptr = 255 : *data_ptr += 127.5;
        }
    }
    cv::imshow("scaling", image_out1);
    cv::waitKey();
    cv::imwrite("scaling.pgm", image_out1);

    return 0;

}
