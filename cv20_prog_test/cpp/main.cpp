#include <iostream>
#include <chrono>
using namespace std;

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

int main() {
    uchar offset = 150;
    double scaling = 2;
    string filename;
    cout << "please input the image name" << endl;

    cin >> filename;
    cv::Mat image = cv::imread("head.pgm", cv::IMREAD_GRAYSCALE);

    bool colorImage;


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
        return -1;
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


//    test.convertTo(test,CV_32F);

    cv::Mat test = image.clone();

//    cv::Mat new_image = cv::Mat::zeros(test.rows,test.cols,test.type());
//    for (size_t y = 0; y < test.rows; y++){
//        for (size_t x = 0; x < test.cols * test.channels(); x++){
////            for(size_t c = 0; c < test.channels(); c++){
//
////                new_image.at<cv::Vec3b>(y,x)[c] =
////                        cv::saturate_cast<uchar>( scaling * test.at<cv::Vec3b>(y,x)[c]);
//
//            auto intensity =  (int)test.at<uchar>(y, x);
////            cout << "\n now " << (int)intensity << endl;
//            intensity -= 127.5 ;
////            cout << " first "<< (int)intensity;
//            intensity *= scaling;
////            cout << " second " << (int)intensity;
//            intensity += 127.5 ;
//            if (intensity < 0) intensity = 0;
//            else if (intensity > 255) intensity = 255;
//            test.at<uchar>(y,x)= (uchar)intensity;
////            }
//        }
//    }

//    visit every pixel with the pointer

    for (size_t y = 0; y < test.rows; y++){
        auto *row_ptr = test.ptr<uchar>(y);
        for(size_t x = 0; x < test.cols * test.channels(); x++) {
            auto *data_ptr = &row_ptr[x];
            auto intensity = (double)*data_ptr;
            intensity -= 127.5;
            intensity *= scaling;
            intensity += 127.5;
            if(intensity > 255) intensity = 255;
            if(intensity < 0) intensity = 0;
            *data_ptr = (uchar)intensity;

        }
    }


    cv::imshow("scaling", test);
    if(!colorImage) cv::imwrite("scaling.pgm", test);
    else cv::imwrite("scaling.ppm", test);
    cv::waitKey();
    cv::destroyAllWindows();

    // automatically called at the end of a scope
    // not needed unless the matrix size varies in different iterations within same loop
    image.release();
    image_out.release();
    test.release();
    return 0;

}
