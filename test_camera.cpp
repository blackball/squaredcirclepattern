#include "squaredcirclepattern.h"
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <string>
#include <vector>

using cv::Point2f;
using std::vector;

// patch the to_string(...) for gcc
namespace patch
{
        template < typename T > std::string to_string( const T& n )
        {
                std::ostringstream stm ;
                stm << n ;
                return stm.str() ;
        }
}

static void drawPoints(cv::Mat &frame, const vector<Point2f> &points)
{
        int face = cv::FONT_HERSHEY_SCRIPT_SIMPLEX;
        
        for (int i = 0; i < points.size(); ++i)
        {
                cv::Point pt(points[i].x, points[i].y);
                cv::circle(frame, pt, 2, cv::Scalar(0,200,0), 2);
                cv::putText(frame, patch::to_string(i), pt, face, 0.6, cv::Scalar(20,20,220), 2);
        }
}

int main(int argc, char *argv[])
{
        SquaredCirclePattern detector;
        cv::VideoCapture cap(0);
        
        if (!cap.isOpened())
                return -1;

        vector<Point2f> points;
        for(;;)
        {
                cv::Mat frame, gray;
                cap >> frame;
                
                if (frame.empty())
                        break;
                cv::cvtColor(frame, gray, cv::COLOR_BGR2GRAY);
                if (0 == detector.detect(gray, points))
                {
                        drawPoints(frame, points);
                }
                cv::imshow("frame", frame);
                if (27 == cv::waitKey(10))
                        break;
        }
        
        return 0;
}
