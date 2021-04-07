#include <iostream>
#include <chrono>
using namespace std;

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

int main() {
    uchar offset = 150;
    float scaling = 2;
    string filename = "baboon.ppm";
//    cout << "please input the image name" << endl;

//    cin >> filename;
//    cv::Mat image = cv::imread("head.pgm", cv::IMREAD_GRAYSCALE);

    bool colorImage;
    cv::Mat image;

    if(filename.substr(filename.find('.') + 1) == "pgm"){
        // read image with only grey scale
        image = cv::imread(filename, cv::IMREAD_GRAYSCALE);
    }else{
        image = cv::imread(filename, cv::IMREAD_COLOR);
    }

    if(image.type() == CV_8UC1){
        colorImage = false;
    }else if(image.type() == CV_8UC3){
        colorImage = true;
    }
//    cout << colorImage<<endl;
    if(image.data == nullptr){
        cerr<<"file doesn't exist"<<endl;
        return 0;
    }

    cout<<"breadth of image is "<< image.cols<<",image height is "<< image.rows <<" number of channel is "
    <<image.channels()<<endl;

    cv::imshow("original image", image);
    cv::waitKey(0);//wait for keyboard input

    //clone the image for output, not directly " = "
    cv::Mat image_out = image.clone();

    //use std::chrono to count the time
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    // visit every pixel with image.at<TYPE>
    for (size_t y = 0; y < image_out.rows; y++){
        for (size_t x = 0; x < image_out.cols * image_out.channels(); x++){
//            for(size_t c = 0; c < image_out.channels(); c++){
                auto intensity =  image_out.at<uchar>(y, x);
//            cout << intensity << endl;
                if (intensity + offset > 255){
                    image_out.at<uchar>(y,x)= 255;
                }
                else{
                    intensity + offset < 0 ?
                    image_out.at<uchar>(y,x) = 0 : image_out.at<uchar>(y,x) += offset;
                }
//            }
        }
    }


    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>> (t2 - t1);
    cout<<"used time is "<< time_used.count() <<endl;
    cv::imshow("output", image_out);
    cv::waitKey();
    if(!colorImage) cv::imwrite("offset.pgm", image_out);
    else cv::imwrite("offset.ppm", image_out);

////     visit every pixel with the pointer
//
//    cv::Mat image_out1 = image.clone();
//    for (size_t y = 0; y < image_out1.rows; y++){
//        auto *row_ptr = image.ptr<uchar>(y);
//        for(size_t x = 0; x < image_out1.cols; x++) {
//            uchar *data_ptr = &row_ptr[x * image_out1.channels()];
//            for(int c = 0; c != image_out1.channels(); c++,data_ptr++){
//                for(int z = 0; z < 10&& x<1 && y <1; z++)
//                    cout << "before:  "<<(int)*data_ptr << endl;
////                cout << "after:   "<< (int)*data_ptr << endl;
//            }
//        }
//    }
//    for (size_t y = 0; y < image_out1.rows; y++){
//        auto *row_ptr = image.ptr<uchar>(y);
//        for(size_t x = 0; x < image_out1.cols; x++) {
//            uchar *data_ptr = &row_ptr[x * image_out1.channels()];
//            for(int c = 0; c != image_out1.channels(); c++,data_ptr++){
////                cout << "before:  "<<(int)*data_ptr << endl;
//                *data_ptr -= 127.5;
//                *data_ptr *= scaling;
//                *data_ptr += 127.5;
//                if(*data_ptr < 0) *data_ptr = 0;
//                if(*data_ptr > 255) *data_ptr = 255;
////                cout << "after:   "<< (int)*data_ptr << endl;
//            }
//        }
//    }
//    ///////////////////////////////////
//    for (size_t y = 0; y < image_out1.rows; y++){
//        auto *row_ptr = image.ptr<uchar>(y);
//        for(size_t x = 0; x < image_out1.cols; x++) {
//            uchar *data_ptr = &row_ptr[x * image_out1.channels()];
//            for(int c = 0; c != image_out1.channels(); c++,data_ptr++){
//                for(int z = 0; z < 10 && x<1 && y <1; z++)
//                    cout << "now:  "<<(int)*data_ptr << endl;
////                cout << "after:   "<< (int)*data_ptr << endl;
//            }
//        }
//    }
//
//    cv::imshow("scaling", image_out1);
//    if(!colorImage) cv::imwrite("scaling.pgm", image_out1);
//    else cv::imwrite("scaling.ppm", image_out1);
//    cv::waitKey();
//    cv::destroyAllWindows();
//
//    // automatically called at the end of a scope
//    // not needed unless the matrix size varies in different iterations within same loop
//    image.release();
//    image_out.release();
//    image_out1.release();
    return 0;

}
