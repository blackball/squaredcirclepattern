/***********************************************************************
 **
 **  Licensed under the Apache License, Version 2.0 (the "License");
 **  you may not use this file except in compliance with the License.
 **  You may obtain a copy of the License at
 **
 **  http://www.apache.org/licenses/LICENSE-2.0
 **
 **  Unless required by applicable law or agreed to in writing, software
 **  distributed under the License is distributed on an "AS IS" BASIS,
 **  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 **  See the License for the specific language governing permissions and
 **  limitations under the License.
 **
 ************************************************************************
 **
 **  Summary :  Squared circle pattern for camera calibration.
 **  Module  :  Squared circle pattern detector.
 **  Author  :  Hui Li
 **  Created :  29-11-2014
 **  Notes   :  Currently, only support one pattern per image.
 **
 ***********************************************************************/


#ifndef SQUAREDCIRCLEPATTERN_H
#define SQUAREDCIRCLEPATTERN_H

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <vector>

class SquaredCirclePatternImpl;

class SquaredCirclePattern
{
public:
        explicit SquaredCirclePattern(const cv::Size sz = cv::Size(8, 5));
	virtual ~SquaredCirclePattern();

        /**
         * Detect pattern points in image file or gray-scale image object.
         * Return 0 if OK, else -1.
         */
        int detect(const char *name, std::vector<cv::Point2f> &out);
        int detect(const cv::Mat &gray, std::vector<cv::Point2f> &out);
        int detect(const IplImage *gray, std::vector<cv::Point2f> &out);

        /**
         * The template points.
         * Return 0.
         */
        int getTemplate(std::vector<cv::Point2f> &out);

        /**
         * Fine tuning utils. No need to use in most cases.
         */
        void setScaling(float scale);
	void setPatternRatio(float minRatio, float maxRatio);
	void setCellRatio(float minRatio, float maxRatio);
	void setMaxCosine(float ratio);
	void setMaxShift(float ratio);
	float getMaxPatternRatio() const;
	float getMinPatternRatio() const;
	float getMaxCellRatio() const;
	float getMinCellRatio() const;
	float getMaxCosine() const;
	float getMaxShift() const;
	float getScaling() const;
private:
	SquaredCirclePatternImpl *m_impl;
};

#endif /* SQUAREDCIRCLEPATTERN_H */
