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

static void test_image(const char *imgName)
{
        cv::Mat im = cv::imread(imgName, 0);
        SquaredCirclePattern detector;
        vector<Point2f> points;
        if (0 == detector.detect(im, points))
        {
                cv::cvtColor(im, im, cv::COLOR_GRAY2BGR);
                drawPoints(im, points);
                cv::imshow("im", im);
                cv::waitKey(0);
        }
        cv::destroyWindow("im");
}

int main(int argc, char *argv[])
{
        if (argc == 2)
        {
                test_image(argv[1]);
        }
        return 0;
}
