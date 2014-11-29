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

#include "squaredcirclepattern.h"
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <vector>
#include <algorithm>
#include <queue>

using std::vector;
using namespace cv;

class SquaredCirclePatternImpl
{
public:
	explicit SquaredCirclePatternImpl(Size patternSize = Size(8,5));
	~SquaredCirclePatternImpl();
	int genTemplate(vector<Point2f> &out) const;
	int detect(const Mat &gray, vector<Point2f> &out);
               
	inline void setPatternRatio(float minRatio, float maxRatio)
        {
                m_maxPatternRatio = maxRatio * m_scale;
                m_minPatternRatio = minRatio * m_scale;
        }
        
	inline void setCellRatio(float minRatio, float maxRatio)
        {
                m_minCellRatio = minRatio;
                m_maxCellRatio = maxRatio;
        }
        
	inline void setMaxCosine(float ratio)
        {
                m_maxCosine = ratio;
        }
        
	inline void setMaxShift(float ratio)
        {
                m_maxShift = ratio;
        }
        
	inline void setScaling(float scale)
        {
                m_scale = scale;
                setPatternRatio(m_minPatternRatio, m_maxPatternRatio);
        }
        
	inline float getMaxPatternRatio() const
        {
                return m_maxPatternRatio;
        }
        
	inline float getMinPatternRatio() const
        {
                return m_minPatternRatio;
        }
        
	inline float getMaxCellRatio() const
        {
                return m_maxCellRatio;
        }
        
	inline float getMinCellRatio() const
        {
                return m_minCellRatio;
        }
        
	inline float getMaxCosine() const
        {
                return m_maxCosine;
        }
        
	inline float getMaxShift() const
        {
                return m_maxShift;
        }
        
	inline float getScaling() const
        {
                return m_scale;
        }

private:
	int findRects(Mat &binary, vector< vector<Point> > &rects, vector<float> &rectAreas);
	vector<Point> findOrder_naive(const Mat &fordraw, vector<Point> &pts, int originI);
private:
	float m_scale;
	float m_minPatternRatio, m_maxPatternRatio;
	float m_minCellRatio, m_maxCellRatio;
	float m_maxCosine, m_maxShift;
	Size m_patternSize;
	vector<Point2f> m_corners;
};

SquaredCirclePattern::SquaredCirclePattern(cv::Size sz)
{
	m_impl = new SquaredCirclePatternImpl(sz);
}

void SquaredCirclePattern::setScaling(float scale)
{
	m_impl->setScaling(scale);
}

void SquaredCirclePattern::setPatternRatio(float minRatio, float maxRatio)
{
	m_impl->setPatternRatio(minRatio, maxRatio);
}

void SquaredCirclePattern::setCellRatio(float minRatio, float maxRatio)
{
	m_impl->setCellRatio(minRatio, maxRatio);
}

void SquaredCirclePattern::setMaxCosine(float mc)
{
	m_impl->setMaxCosine(mc);
}

void SquaredCirclePattern::setMaxShift(float ms)
{
	m_impl->setMaxShift(ms);
}

float SquaredCirclePattern::getMaxCellRatio() const
{
	return m_impl->getMaxCellRatio();
}

float SquaredCirclePattern::getMinCellRatio() const
{
	return m_impl->getMinCellRatio();
}

float SquaredCirclePattern::getMaxPatternRatio() const
{
	return m_impl->getMaxPatternRatio();
}

float SquaredCirclePattern::getMinPatternRatio() const
{
	return m_impl->getMinPatternRatio();
}

float SquaredCirclePattern::getMaxCosine() const
{
	return m_impl->getMaxCosine();
}

float SquaredCirclePattern::getMaxShift() const
{
	return m_impl->getMaxShift();
}

float SquaredCirclePattern::getScaling() const
{
	return m_impl->getScaling();
}

SquaredCirclePattern::~SquaredCirclePattern()
{
	delete m_impl;
}

int SquaredCirclePattern::getTemplate(vector<Point2f> &out)
{
        return m_impl->genTemplate(out);
}

int SquaredCirclePattern::detect(const char *name, vector<Point2f> &out)
{
	Mat gray = imread(name, 0); // gray-scale 
	if (gray.empty())
	{
		return -1;
	}
	return m_impl->detect(gray, out);
}

int SquaredCirclePattern::detect(const IplImage *gray, vector<Point2f> &out)
{
        cv::Mat m = cv::Mat(gray, true);
        return m_impl->detect(m, out);
}

int SquaredCirclePattern::detect(const Mat &gray, vector<Point2f> &out)
{
	return m_impl->detect(gray, out);
}

/****Implementation*****/
SquaredCirclePatternImpl::SquaredCirclePatternImpl(Size patternSize)
{
	m_patternSize = patternSize;
	// generate corners
	const float X = (patternSize.width - 1) * 100.f;
	const float Y = (patternSize.height - 1) * 100.f;
	m_corners.push_back(Point2f(0, 0));
	m_corners.push_back(Point2f(0, Y));
	m_corners.push_back(Point2f(X, 0));
	m_corners.push_back(Point2f(X, Y));

	m_scale = 1.0f;
	m_minPatternRatio = 0.01f * m_scale;
	m_maxPatternRatio = 0.8f * m_scale;
	m_maxCosine = 0.45;
	//m_maxShift = 0.45;
	m_minCellRatio = 0.05;
	m_maxCellRatio = 2.0;
}

SquaredCirclePatternImpl::~SquaredCirclePatternImpl()
{
	
}

int SquaredCirclePatternImpl::genTemplate(vector<Point2f> &points) const
{
        points.resize(m_patternSize.height * m_patternSize.width);
	for (int i = 0; i < m_patternSize.height; ++i)
	{
		for (int j = 0; j < m_patternSize.width; ++j)
		{
			points[i*m_patternSize.width + j] = Point2f(j*100.f, i*100.f);
		}
	}
	return 0;
}

static inline double ellipse_area(const RotatedRect &rr)
{
	return rr.size.width * rr.size.height * 0.25 * 3.141593;
}

static inline double angle(Point pt1, Point pt2, Point pt0 ) 
{
        const double dx1 = pt1.x - pt0.x;
        const double dy1 = pt1.y - pt0.y;
        const double dx2 = pt2.x - pt0.x;
        const double dy2 = pt2.y - pt0.y;
        return (dx1*dx2 + dy1*dy2)/sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10);
}

static inline double distanceSquare(const Point a, const Point b)
{
	return 1.0 * ((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

static inline double distanceP(const Point a, const Point b)
{
	return sqrt((double)(a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

/* @TODO The better way is to using tracing algorithm, this naive algorithm is not robust to noisy image */
static vector<Point> findEdgePoints(const Mat &im, Point start, Rect boundingRect)
{
	// 0. find the starting edge point
	const int x_start = boundingRect.x;
	const int y_start = boundingRect.y;
	const int x_end = boundingRect.x + boundingRect.width;
	const int y_end = boundingRect.y + boundingRect.height;
	vector<Point> edgePoints;
	int y = start.y, x;
	while (y > y_start)
	{
		Point left, right;		
		// left point
		for (x = start.x; x > x_start; --x)
		{
			if (im.at<uchar>(y,x) == 0)
			{
				left.x = x;
				left.y = y;
				break;
			}
		}
                
		assert(x != x_start);
		// right point 
		for (x = start.x; x < x_end; ++x)
		{
			if (im.at<uchar>(y,x) == 0)
			{
				right.x = x;
				right.y = y;
				break;
			}
		}
                
		assert(x != x_end);
		if (left.x >= right.x)
		{
			break;
		}
		edgePoints.push_back(left);
		edgePoints.push_back(right);
		--y;
	}
	
	y = start.y+1;
	while (y < y_end)
	{
		Point left, right;

                // left point
		for (x = start.x; x > x_start; --x)
		{
			if (im.at<uchar>(y,x) == 0)
			{
				left.x = x;
				left.y = y;
				break;
			}
		}
                
		assert(x != x_start);

                // right point 
		for (x = start.x; x < x_end; ++x)
		{
			if (im.at<uchar>(y,x) == 0)
			{
				right.x = x;
				right.y = y;
				break;
			}
		}
                
		assert(x != x_end);
                
		if (left.x >= right.x)
		{
			break;
		}
                
		edgePoints.push_back(left);
		edgePoints.push_back(right);
		++y;
	}
	return edgePoints;
}

static inline Point meanPoint(const vector<Point> &pts)
{
	Point mp(0,0);
	for (size_t i = 0; i < pts.size(); ++i) 
		mp += pts[i];
	mp.x /= pts.size();
	mp.y /= pts.size();
	return mp;
}

static inline void removePoint(vector<Point> &corners, int i0)
{
	corners[i0] = corners[corners.size() - 1];
	corners.pop_back();
}

static inline int getNNOrder(const Point2f p, int &x, int &y)
{
	x = cvRound(p.x/100.0);
	y = cvRound(p.y/100.0);
	if (x >= 0 && x <= 7 && y >= 0 && y <= 4)
	{
		return 0;
	}
	return -1;
}

static void findCorners(const vector<Point> &pts, vector<Point> &corners)
{
	vector<Point> hull;
	convexHull(pts, hull, true);
        
	const int hsz = hull.size();
	// find four ~90 corners
	vector<double> angles;
        
	angles.push_back(-1.0 * fabs(angle(hull[hsz - 1], hull[1], hull[0])));
	for (size_t i = 1; i < hull.size() - 1; ++i)
	{
		angles.push_back(-1.0 * fabs(angle(hull[i-1], hull[i+1], hull[i])));
	}
	angles.push_back(-1.0 * fabs(angle(hull[hsz-2], hull[0], hull[hsz-1])));

	// naive way to find smallest 4 points
	std::priority_queue<std::pair<double, int> > pq;
	for (size_t i = 0; i < angles.size(); ++i)
	{
		pq.push(std::pair<double, int>(angles[i], i));
	}
	
	for (int i = 0; i < 4; ++i)
	{
		corners[i] = hull[ pq.top().second ];
		pq.pop();
	}
}

static inline int nearestPoint(const vector<Point> &pts, Point pt)
{
	int ni = 0;
	double nd = distanceSquare(pt, pts[0]);
	for (size_t i = 1; i < pts.size(); ++i)
	{
		const double td = distanceSquare(pt, pts[i]);
		if (td < nd)
		{
			nd = td;
			ni = i;
		}
	}
	return ni;
}

vector<Point> SquaredCirclePatternImpl::findOrder_naive(const Mat &binary, vector<Point> &pts, int originI)
{
	// a naive approach 
	// 0 find four corners
	// 1 map the corners
	// 2 transform points
	// 3 ordering
	vector<Point> corners(4);
	findCorners(pts, corners);
	Point op = pts[originI];
	
	// map the corners
	Point p0,p1,p2,p3;
	int i0 = nearestPoint(corners, op);
	p0 = corners[i0];
	removePoint(corners, i0);
	int i1 = nearestPoint(corners, p0);
	p1 = corners[i1];
	removePoint(corners, i1);
	int i2 = nearestPoint(corners, p0);
	int i2_1 = nearestPoint(corners, p1);
	if (i2 != i2_1) {
		p2 = corners[i2];
		removePoint(corners, i2);
		p3 = corners[0];
	}
	else {
		p2 = corners[1 - i2];
		removePoint(corners, 1 - i2);
		p3 = corners[0];
	}

	corners.resize(4);
	corners[0] = p0;
	corners[1] = p1;
	corners[2] = p2;
	corners[3] = p3;

	vector<Point2f> cornersf(4);
        for (int i = 0; i < cornersf.size(); ++i)
        {
                cornersf[i].x = corners[i].x;
                cornersf[i].y = corners[i].y;
        }
        
	vector<Point2f> &tplCorners = m_corners;
	vector<Point2f> ptsT(pts.size());

	// transform those points 
	Mat H = getPerspectiveTransform(cornersf, tplCorners);
	vector<Point2f> ptsf(pts.size());
        for (int i = 0; i < ptsf.size(); ++i)
        {
                ptsf[i].x = pts[i].x;
                ptsf[i].y = pts[i].y;
        }
        
        perspectiveTransform(ptsf, ptsT, H);

	vector<Point> orderedPts(pts.size());
	// find the right order 
	// For point in ptsT, find the neaset points in template, then 
	// save order
	for (size_t i = 0; i < ptsT.size(); ++i)
	{
		int x, y;
		if (0 == getNNOrder(ptsT[i], x, y))
                {
			orderedPts[x+y*8] = pts[i];
		}
		else
                {
			return vector<Point>();
		}
	}
	return orderedPts;
}

int SquaredCirclePatternImpl::detect(const Mat &gray, vector<Point2f> &out)
{
	Mat binary, smallBinary;
	threshold(gray, binary, 0, 255, CV_THRESH_OTSU | CV_THRESH_BINARY);
        
        if (fabs(m_scale - 1.0) < 0.1)
        {
                smallBinary = binary.clone();
        }
        else
        {
                resize(binary, smallBinary, Size(binary.cols * m_scale, binary.rows * m_scale));
        }
	
	vector< vector<Point> > rects;
	vector<float> rectAreas;

	if (-1 == findRects(smallBinary, rects, rectAreas))
	{
		return -1;	
	}

	if (rects.size() != m_patternSize.width * m_patternSize.height)
	{
		return -1;
	}

	float inv_scale = 1.0 / m_scale;
	// scale the rects and areas back
	for (size_t i = 0; i < rects.size(); ++i)
	{
		for (size_t j = 0; j < rects[i].size(); ++j)
		{
			rects[i][j] *= inv_scale;
		}
	}
	for (size_t i = 0; i < rectAreas.size(); ++i)
	{
		rectAreas[i] *= inv_scale * inv_scale;
	}

	vector<Point> centers;
	double maxRatio = 0;
	Point originP;
	int originI = 0;
	for (size_t i = 0; i < rects.size(); ++i)
	{
		Rect br = boundingRect(rects[i]);
		Point start = meanPoint(rects[i]);
		vector<Point> edgePoints = findEdgePoints(binary, start, br);
		if (edgePoints.size()) {
			RotatedRect rr = fitEllipse(edgePoints);
			// calculate ellipse area
			double rea = ellipse_area(rr) / rectAreas[i];
			if (rea > maxRatio)
			{
				maxRatio = rea;
				originP = rr.center;
				originI = i;
			}
			centers.push_back(rr.center);
		}
	}

	vector<Point> orderedPts = findOrder_naive(binary, centers, originI);

	if (orderedPts.size() > 0)
        {
		out.resize(orderedPts.size());
                for (int i = 0; i < out.size(); ++i)
                {
                        out[i].x = orderedPts[i].x;
                        out[i].y = orderedPts[i].y;
                }
		return 0;
	}
	return -1;
}

static inline vector< vector<int> > sumInternalContours(const vector<Vec4i> &H)
{
	vector< vector<int> > v(H.size(), vector<int>());
	for (size_t i = 0; i < H.size(); ++i)
	{
		int t = H[i][3];
		while (t != -1)
		{
			v[t].push_back(i);
			t = H[t][3];
		}
	}
	return v;
}

int SquaredCirclePatternImpl::findRects(Mat& bim, vector< vector<Point> >& rects, vector<float>& rectAreas)
{
	vector< vector<Point> > contours;
	vector<Point> approx;
	Mat dbim;
        
	// Use median filter to remove some noise 
	medianBlur(bim, dbim, 1);

	vector<Vec4i> hierarchy;
	findContours(dbim, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE);
	vector< vector<int> > internalContours = sumInternalContours(hierarchy);

	const float imarea = bim.size().area();

	int targetI = -1;
	int maxKidNum = 0; // Here maxKidNum is not the meaning of its name indicates.
	int exactNum = this->m_patternSize.height * this->m_patternSize.width * 2;
	float parentArea = 0;
	
	for (size_t i = 0; i < contours.size(); ++i)
	{
		int kidNum = internalContours[i].size();
		if ( kidNum < m_patternSize.width * m_patternSize.height)
		{
			continue;
		}
		else
		{
			approxPolyDP(Mat(contours[i]), approx, arcLength(Mat(contours[i]), true) * 0.02, true);
			
                        if (approx.size() >= 4 && approx.size() <= 8)
			{
				RotatedRect rrect = minAreaRect(contours[i]);
				const float area = rrect.size.area();
				if ( area > imarea * m_minPatternRatio && area < imarea * m_maxPatternRatio)
				{
					if (kidNum > (exactNum/2))
					{
						int sa = abs(kidNum - exactNum);
						int sb = abs(maxKidNum - exactNum);
						if (sa < sb) {
							maxKidNum = kidNum;
							targetI = i;
							parentArea = area;
						}
					}
				}
			}
		}
	}

	if (targetI == -1)
	{
		return -1;
	}

	for (size_t i = 0; i < internalContours[targetI].size(); ++i)
	{
		const int id = internalContours[targetI][i];
		approxPolyDP(Mat(contours[id]), approx, arcLength(Mat(contours[id]), true) * 0.02, true);
		const cv::Mat mm = Mat(approx);
		if (approx.size() == 4 && isContourConvex(Mat(approx)))
		{
			RotatedRect rrect = minAreaRect(contours[id]);
			const float area = rrect.size.area();//fabs(contourArea(mm));
			const float minAreaBound = parentArea * m_minCellRatio / (m_patternSize.width * m_patternSize.height);
			const float maxAreaBound = parentArea * m_maxCellRatio / (m_patternSize.width * m_patternSize.height);
			if (area >  minAreaBound && area < maxAreaBound)
			{
				double maxCosine = 0.0;
                                // angle shift was disabled here for it's in-effectiveness.
				//double maxShift = 0;
				for (int j = 2; j < approx.size() + 1; j++)
				{
					double cosine = fabs(angle(approx[j%4], approx[j-2], approx[j-1]));
					maxCosine = MAX(maxCosine, cosine);
					//double da = distanceP(approx[j%4], approx[j-1]);
					//double db = distanceP(approx[j-2], approx[j-1]);
					//da = fabs(da/db - 1.0);
					//maxShift = MAX(maxShift, da);
				}

				if (maxCosine < m_maxCosine /* && maxShift < m_maxShift*/) 
				{
					rects.push_back(approx);
					rectAreas.push_back(area);
				}
			}
		}	
	}

	if (rects.size() != m_patternSize.width * m_patternSize.height)
	{
		return -1;
	}

	return 0;
}
