//============================================================================
// Name        : sample2d.cpp
// Author      : Stefan Ohrhallinger
// Version     :
// Copyright   : GPL v3 or later
// Description : Sample a cubic Bezier curve with epsilon-sampling
//============================================================================

#define QUINTICSOLVER
#define NANOSVG_IMPLEMENTATION

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include <cassert>
#include <unsupported/Eigen/Polynomials>
#include <ANN/ANN.h>
#include "PrecTimer.h"
#include "nanosvg.h"
#include "root_finder.hpp"

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>

using namespace std;
using namespace casa;

#define PI 3.1415926
#define SQR(a) ((a)*(a))
template<typename T> inline T CUB(T t) { return (SQR(t)*(t)); };

const double EPSILON_T = 0.000000001;	// precision of point on curve

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point CGPoint;
typedef CGAL::Triangulation_2<K, Tds> Triangulation_2;

struct Point
{
	array<float, 2> data;

	Point() {}

	Point(const float& x, const float& y)
	{
		data[0] = x;
		data[1] = y;
	}

	float & operator[](int ind) {return data[ind];}

	Point operator-(const Point &rhs) const
	{
		return Point(this->x()-rhs.x(), this->y()-rhs.y());
	}

	Point operator+(const Point &rhs) const
	{
		return Point(this->x()+rhs.x(), this->y()+rhs.y());
	}

	//dot product
	float operator*(const Point &rhs) const
	{
		return (this->x()*rhs.x()) + (this->y()*rhs.y());
	}

	//scale
	Point operator*(const float &s) const
	{
		return Point(this->x()*s, this->y()*s);
	}

	// dot product
	float dot(Point p)
	{
		return x()*p.x() + y()*p.y();
	}

	float squared_length()
	{
		return x()*x() + y()*y();
	}

	float squared_distance(const Point &p0)
	{
		auto d = p0 - *this;
		return d.squared_length();
	}

	float distance(const Point &p0)
	{
		return sqrt(squared_distance(p0));
	}

	void normalize()
	{
		float l = sqrt(squared_length());
		data[0] /= l;
		data[1] /= l;
	}

	const float &x() const  {return data[0];}
	const float &y() const  {return data[1];}
};

typedef struct
{
	double x, y;
} CBPoint;

typedef struct
{
	int v[2];
} CBEdge;

struct DPoint
{
	array<double, 2> data;

	DPoint() {}

	DPoint(const double& x, const double& y)
	{
		data[0] = x;
		data[1] = y;
	}

	DPoint(const CBPoint p)
	{
		data[0] = p.x;
		data[1] = p.y;
	}

	double & operator[](int ind) { return data[ind]; }

	DPoint operator-(const DPoint &rhs) const
	{
		return DPoint(this->x() - rhs.x(), this->y() - rhs.y());
	}

	DPoint operator+(const DPoint &rhs) const
	{
		return DPoint(this->x() + rhs.x(), this->y() + rhs.y());
	}

	// dot product
	double operator*(const DPoint &rhs) const
	{
		return (this->x()*rhs.x()) + (this->y()*rhs.y());
	}

	//scale
	DPoint operator*(const double &s) const
	{
		return DPoint(this->x()*s, this->y()*s);
	}

	// dot product
	double dot(DPoint p)
	{
		return x()*p.x() + y()*p.y();
	}

	double squared_length()
	{
		return x()*x() + y()*y();
	}

	double length()
	{
		return sqrt(squared_length());
	}

	double squared_distance(const DPoint &p0)
	{
		auto d = p0 - *this;
		return d.squared_length();
	}

	double distance(const DPoint &p0)
	{
		return sqrt(squared_distance(p0));
	}

	void normalize()
	{
		double l = sqrt(squared_length());
		data[0] /= l;
		data[1] /= l;
	}

	const double &x() const { return data[0]; }
	const double &y() const { return data[1]; }
};

/*
 * return distance of p0 from line segment p1-p2
 * project p1-p0 onto p1-p2 (normalized) to find t, then compute n on p1-p2 or closest end point if out of range
 */
double distancePFromLineSegment(DPoint p0, DPoint p1, DPoint p2)
{
	DPoint v = p2 - p1;
	double len = v.length();
	double t = ((p0[0] - p1[0])*v[0] + (p0[1] - p1[1])*v[1])/SQR(len);

	if ((t >= 0.0) && (t <= 1.0))
	{
		DPoint n = p1 + v*t;

		return n.distance(p0);
	}
	else
	{
		double dist[2] = { p0.squared_distance(p1), p0.squared_distance(p2) };

		return (dist[0] < dist[1]) ? sqrt(dist[0]) : sqrt(dist[1]);
	}
}

/*
 * evaluate cubic Bezier curve for parameter t
 */
DPoint evalCBC(DPoint *cbc, double t)
{
	return cbc[0]*CUB(1 - t) + cbc[1]*3*t*SQR(1 - t) + cbc[2]*3*SQR(t)*(1 - t) + cbc[3]*CUB(t);
}

/*
 * determine closest point in cubic Bezier curve to point p
 */
DPoint closestPointOnCBC(DPoint *cbc, DPoint p, double &t)
{
	int i, j, k;

	// translate curve points to p as origin
	DPoint b[4];

	for (i = 0; i < 4; i++)
		b[i] = cbc[i] - p;

	// f(t) = |B3(t)-q|^2 where q=(0,0) transformed to origin, so |B3(t)|^2
	// set f'(t) = 0
	Eigen::MatrixXd bezier(4, 4);	// 4x4 Bezier matrix M4

	bezier << -1, 3, -3, 1,
			  3, -6, 3, 0,
			  -3, 3, 0, 0,
			  1, 0, 0, 0;

	// multiply Bezier matrix with P to get B3 (both coordinates)
	Eigen::MatrixXd coeffs42 = Eigen::MatrixXd::Zero(4, 2);	// B4=M4*P

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 2; k++)
				coeffs42(i, k) += bezier(i, j)*b[j][k];

	// then multiply polynomial vectors for each coordinate
	Eigen::MatrixXd coeffs72 = Eigen::MatrixXd::Zero(7, 2);	// B4^2

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 2; k++)
				coeffs72(6 - (i + j), k) += coeffs42(3 - i, k)*coeffs42(3 - j, k);

	Eigen::VectorXd coeffs7 = Eigen::VectorXd::Zero(7);	// sum up coordinates to coefficients of polynomial

	for (i = 0; i < 7; i++)
		for (k = 0; k < 2; k++)
			coeffs7(i) += coeffs72(i, k);

	Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(6);	// compute (B4^2)'

#ifdef QUINTICSOLVER
	for (i = 0; i < 6; i++)
		coeffs(i) = (6 - i)*coeffs7(i);	// solver needs highest coefficient ordered first

	double QUINTICSOLVER_TOLERANCE = 0.000000001;
	set<double> roots = RootFinder::solvePolynomial(coeffs, 0.0, 1.0, QUINTICSOLVER_TOLERANCE);
#else
	for (i = 0; i < 6; i++)
		coeffs(5 - i) = (6 - i)*coeffs7(i);	// solver needs highest coefficient ordered last

	Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
	solver.compute(coeffs);
	list<double> roots;
	solver.realRoots(roots);
#endif

	// determine closest point for real roots: if no roots, check which end point is closest, else test for all roots
	DPoint q;

	double sqDist = numeric_limits<double>::max();

	for (auto r:roots)
		if ((r >= 0.0) && (r <= 1.0))
		{
			DPoint qq = evalCBC(b, r);

			if (qq.squared_length() < sqDist)
			{
				sqDist = qq.squared_length();
				q = qq;
				t = r;
			}
		}

	int closerIndex = (b[0].squared_length() < b[3].squared_length()) ? 0 : 1;

	if ((roots.size() == 0) || (b[closerIndex*3].squared_length() < sqDist))
	{
		q = b[closerIndex*3];
		t = closerIndex;
	}

	return p + q;	// add point back
}

/*
 * return radius for circle through point p with normalized normal n and point q
 * q can be mirrored on the line through n, therefore the radius is the circumradius of the triangle pqq'
 */
double radiusForCircleThrough2PointsandNormal(DPoint p, DPoint n, DPoint q)
{
	double a2 = p.squared_distance(q);
	double a = sqrt(a2);
	DPoint n2(-n[1], n[0]);
	DPoint pq = p - q;
	double dist = abs(pq*n2);
	double b = 2*dist;

	if (b == 0)
		return 0.5*sqrt(pq.squared_length());	// distance pq = diameter

	double d2 = (2*a + b)*SQR(b)*(2*a - b);	// r=abc/sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) -> isosceles a=b

	if (d2 <= 0.0)
		return numeric_limits<double>::max();	// triangle points are collinear, infinite radius

	double d = sqrt(d2);
	double r = abs(a2*b/d);

	return r;	// circumradius of triangle adapted to isosceles version
}

int iterateMedialPoint(vector<vector<DPoint> > &cbcs, DPoint currP, DPoint normal, double &minR, int sampleCount)
{
	const double DIST_EPSILON = 0.000000001;
	int i, j;
	DPoint mPoint = currP + normal*minR;

	// successively compute more accurate medial point by determining closest point of all curve segments
	bool adjusted;
	double minDist;

	// DEBUG
	int iterations = 0;

	do
	{
		adjusted = false;
		minDist = minR;
		DPoint minN;

		for (auto cbc:cbcs)
		{
			for (i = 0; i < (int)cbc.size()/4; i++)
			{
				// test with convex hull of bezier curve first:

				// first if outside
				DPoint centroid(0.0, 0.0), normal[4];

				for (j = 0; j < 4; j++)
				{
					centroid = centroid + cbc[i*4 + j];
					DPoint edge = cbc[i*4 + ((j + 1) % 4)] - cbc[i*4 + j];
					normal[j][0] = edge[1];
					normal[j][1] = -edge[0];
				}

				centroid = centroid*0.25;
				bool insideOrient = normal[0].dot(centroid - cbc[i*4]) < 0.0;
				bool isInside = true;
				j = 0;

				while (isInside && (j < 4))
				{
					isInside = (normal[j].dot(mPoint - cbc[i*4 + j]) < 0.0) == insideOrient;
					j++;
				}

				bool isHullCloser = isInside;

				// all edges between vertices in case vertices are not ordered in convex hull
				int edges[6][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 0, 2 }, { 1, 3 } };

				// then if convex hull is closer to current mPoint
				j = 0;
				while (!isHullCloser && (j < 6))
				{
					isHullCloser = (distancePFromLineSegment(mPoint, cbc[i*4 + edges[j][0]], cbc[i*4 + edges[j][1]]) - DIST_EPSILON < minDist);
					j++;
				}

				if (isInside || isHullCloser)
				{
					double curveT;
					DPoint n = closestPointOnCBC(&cbc[i*4], mPoint, curveT);
					double dist = n.distance(mPoint);

					if ((dist < minDist) && fmod((i + curveT)*sampleCount, (cbc.size()/4)*sampleCount) >= 1.0)
					{
						minN = n;
						minDist = dist;
					}
				}
			}
		}

		if (minDist + DIST_EPSILON < minR)	// iterate to desired precision
		{
			minR = radiusForCircleThrough2PointsandNormal(currP, normal, minN);
			mPoint = currP + normal*minR;
			adjusted = true;
		}

		// DEBUG
		iterations++;

	} while (adjusted);

	return iterations;
}

/*
 * compute distance to piece-wise linear medial axis
 */
double computeDistToMedialAxis(DPoint p, vector<vector<DPoint> > &cbcs, vector<DPoint> &curvePoints, vector<DPoint> &mPoints, vector<bool> &orient, int curveIndex, ANNkd_tree *kdtree, int sampleCount, DPoint prevP)
{
	// first, get nearest neighbor in points sampled on medial axis
	ANNpointArray search_point = annAllocPts(1, 2);
	ANNidxArray nnIdx = new ANNidx[1];
	ANNdistArray distances = new ANNdist[1];

	search_point[0][0] = p[0];
	search_point[0][1] = p[1];
	kdtree->annkSearch(search_point[0], 1, nnIdx, distances);
	double distP = sqrt(distances[0]);

	double distM = distP;

	// then, test whether incident piece-wise linear segments to medial point intersect circle with radius distP
	int curveSize = cbcs[curveIndex].size()/4*sampleCount;
	int currIndex = nnIdx[0];
	int prevIndex = (currIndex + curveSize - 1) % curveSize;
	int nextIndex = (currIndex + 1) % curveSize;
	DPoint currM(mPoints[currIndex]);

	if ((&orient != NULL) && (orient.size() > 0))
	{
		if (orient[currIndex] == orient[prevIndex])
		{
			double dist = distancePFromLineSegment(p, mPoints[prevIndex], currM);

			if (abs(dist - distM) < p.distance(prevP))	// sanity check using 1-Lipschitz
				if (dist < distP)
					distP = dist;
		}

		if (orient[currIndex] == orient[nextIndex])
		{
			double dist = distancePFromLineSegment(p, currM, mPoints[nextIndex]);

			if (abs(dist - distM) < p.distance(prevP))	// sanity check using 1-Lipschitz
				if (dist < distP)
					distP = dist;
		}
	}

	return distP;
}

/*
 * compute LFS values at specified curve point
 */
double computeLFSForCurveAt(vector<vector<DPoint> > &cbcs, vector<DPoint> &curvePoints, vector<DPoint> &mPoints, vector<bool> &orient, int curveIndex, double t, double prevT, ANNkd_tree *kdtree, int sampleCount)
{
	// compute point on curve
	DPoint prevP = evalCBC(&cbcs[curveIndex][4*(int)(prevT/sampleCount)], fmod(prevT, sampleCount)/sampleCount);
	DPoint p = evalCBC(&cbcs[curveIndex][4*(int)(t/sampleCount)], fmod(t, sampleCount)/sampleCount);
	double distP = computeDistToMedialAxis(p, cbcs, curvePoints, mPoints, orient, curveIndex, kdtree, sampleCount, prevP);

	return distP;
}

/*
 * compute normal for curve
 */
DPoint computeNormalForCurve(double t, DPoint *cb)
{
	int i;
	Point qb[3];
	DPoint tangent;

	// compute tangent as derivative

	for (i = 0; i < 2; i++)
	{
		qb[0][i] = 3*(cb[1][i] - cb[0][i]);
		qb[1][i] = 3*(cb[2][i] - cb[1][i]);
		qb[2][i] = 3*(cb[3][i] - cb[2][i]);
	}

	for (i = 0; i < 2; i++)
		tangent[i] = qb[0][i]*SQR(1 - t) + 2*qb[1][i]*t*(1 - t) + qb[2][i]*SQR(t);

/*
	// -3(1-t)^2 * P0 + 3(1-t)^2 * P1 - 6t(1-t) * P1 - 3t^2 * P2 + 6t(1-t) * P2 + 3t^2 * P3
	for (i = 0; i < 2; i++)
		tangent[i] = -3*SQR(1 - t)*cb[0][i] + (3*SQR(1 - t) - 6*t*(1 - t))*cb[1][i] + (-3*SQR(t) + 6*t*(1 - t))*cb[2][i] + 3*SQR(t)*cb[3][i];
*/

	// handle quadratic Bezier (cb[0] == cb[1])
	if ((tangent[0] == 0.0) && (tangent[1] == 0.0))
	{
		tangent[0] = cb[2][0] - cb[0][0];
		tangent[1] = cb[2][1] - cb[0][1];
	}

	// normalize
	tangent.normalize();

	// hack for very close points
	if (isnan(tangent[0]) || isnan(tangent[1]))
	{
		tangent[0] = 1.0;
		tangent[1] = 0.0;
	}

	// then rotate by 90 degrees
	return DPoint(tangent[1], -tangent[0]);
}

/*
 * compute medial points for all curve points
 */
void computeMedialPoints(vector<vector<DPoint> > &cbcs, vector<DPoint> &curvePoints, vector<bool> &orient, int curveIndex, vector<DPoint> &medialPoints, int sampleCount)
{
	int i, j;
	DPoint mPoint;
	bool newOrient = false;

	// compute bounding box extent
	DPoint min(numeric_limits<double>::max(), numeric_limits<double>::max());
	DPoint max(-numeric_limits<double>::max(), -numeric_limits<double>::max());

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		DPoint m = curvePoints[i];

		for (j = 0; j < 2; j++)
		{
			if (m[j] < min[j])
				min[j] = m[j];

			if (m[j] > max[j])
				max[j] = m[j];
		}
	}

	double bbDiag = min.distance(max);
	orient.resize(curvePoints.size());

	for (i = 0; i < (int)curvePoints.size(); i++)
	{
		curvePoints[i] = evalCBC(&cbcs[curveIndex][4*(int)(i/sampleCount)], (double)(i % sampleCount)/sampleCount);
		DPoint currP = curvePoints[i];

		DPoint normal = computeNormalForCurve((double)(i % sampleCount)/sampleCount, &cbcs[curveIndex][4*(int)(i/sampleCount)]);
		DPoint oppNormal = normal;
		oppNormal[0] = -normal[0];
		oppNormal[1] = -normal[1];
		double minR[2] = { 0.5*bbDiag, 0.5*bbDiag };

		// approximate radius from next discrete point on same curve, valid for a single side
		DPoint curr2P = curvePoints[(i + 1) % curvePoints.size()];
		double minRNext = radiusForCircleThrough2PointsandNormal(currP, normal, curr2P);
		double dot = normal.dot(curr2P - currP);	// determine side for that radius
		minR[(dot > 0.0) ? 0 : 1] = minRNext;

		int iterations = 0;

		if (minR[0] != numeric_limits<double>::max())
			iterations += iterateMedialPoint(cbcs, currP, normal, minR[0], sampleCount);

		if (minR[1] != numeric_limits<double>::max())
			iterations += iterateMedialPoint(cbcs, currP, oppNormal, minR[1], sampleCount);

		newOrient = (minR[1] < minR[0]);

		if (newOrient)
			mPoint = currP - normal*minR[1];
		else
			mPoint = currP + normal*minR[0];

		medialPoints.push_back(mPoint);

		orient[i] = newOrient;
	}
}

/*
 * return distance of p0 from line p1-p2
 */
double distancePFromLine(DPoint p0, DPoint p1, DPoint p2)
{
	DPoint normal(p2 - p1);
	swap(normal[0], normal[1]);
	normal[0] = -normal[0];
	normal.normalize();

	return abs(normal*(p0 - p1));
}

/*
 * generate point set (+edges) by sampling with condition < epsMax, max error to original curve and optionally manifold condition
 */
void sampleCurvesByEps(vector<pair<vector<DPoint>, bool> > &curves, vector<vector<DPoint> > &cbcs, float maxEps,
		float error, float perturb, vector<vector<CBPoint> > &samples, vector<vector<double> > &sValues,
		vector<vector<float> > &noise, float offset, vector<DPoint> &medialPoints, vector<vector<bool> > &orient,
		vector<vector<double> > &lfsCurve, int sampleCount, float minDist, float maxDist)
{
	int i, j, k, pointCount = 0;

	// TEST
	vector<double> s;

	// for each connected curve, count points
	for (auto curve:curves)
		pointCount += curve.first.size();

	ANNpointArray ann_points;
	ann_points = annAllocPts(pointCount, 2);

	// for each connected curve, compute their medial points
	int index = 0;

	for (k = 0; k < (int)curves.size(); k++)
	{
		vector<DPoint> *curve = &curves[k].first;
		computeMedialPoints(cbcs, *curve, orient[k], k, medialPoints, sampleCount);

		for (i = 0; i < (int)curve->size(); i++)
		{
			auto q = ann_points[index];
			q[0] = medialPoints[index][0];
			q[1] = medialPoints[index][1];
			index++;
		}

		orient[k].resize(curve->size());
	}

	ANNkd_tree *kdtree = new ANNkd_tree(ann_points, pointCount, 2);

	// for each connected curve, sample with epsilon
	for (k = 0; k < (int)curves.size(); k++)
	{
		vector<DPoint> *curve = &curves[k].first;
		bool isClosed = curves[k].second;
		lfsCurve[k].resize(curve->size());

		// compute lfs of the curve points
		for (i = 0; i < (int)curve->size(); i++)
			lfsCurve[k][i] = computeLFSForCurveAt(cbcs, *curve, medialPoints, orient[k], k, i, (i + curve->size() - 1) % curve->size(), kdtree, sampleCount);

		vector<Point> offsetCurve(curve->size());
		int intOfs = (int)(curve->size()*offset) % curve->size();
		i = intOfs;
		double prevT = 0.0, t = intOfs;
		int prevI, prevJ = 0;
		DPoint prevSample, candSample, newSample = (*curve)[i];
		bool isCompleted = false;

		do
		{
			// add new sample
			DPoint normal = computeNormalForCurve(fmod(t, sampleCount)/sampleCount, &cbcs[k][4*((int)t/sampleCount)]);

			if (perturb != 0.0)
			{
				// perturb new sample by up to perturb*lfs(p) along its normal in any direction
				double random = (double)rand()/RAND_MAX*2 - 1.0;
				newSample = newSample + normal*(random*perturb*lfsCurve[k][round(t)]);
			}

			CBPoint cbSample;
			cbSample.x = newSample[0];
			cbSample.y = newSample[1];
			samples[k].push_back(cbSample);
			sValues[k].push_back(t);
			noise[k].push_back(perturb*lfsCurve[k][round(t)]);

			// TEST
			s.push_back(t);

			prevSample = newSample;
			prevT = t;

			// move to next point
			prevI = i;
			i = floor(t) + 1;

			if (i == (int)curve->size())
				i = 0;

			// test candidate sample if it conforms to the sampling condition
			bool isConforming = true;

			// increase i to next point x while epsilon*lfs(x)<|prevSample,x| holds
			while (!isCompleted && isConforming)
			{
				if (i == intOfs)
					isCompleted = true;

				DPoint candX = (*curve)[i];
				isConforming = true;

				if (isConforming)
				{
					// test for epsilon condition (a sample within dist/lfs < maxEps)
					double lfs = lfsCurve[k][i];
					double dist = candX.distance(prevSample);
					isConforming = dist/lfs < maxEps;
				}

				if (isConforming)
				{
					prevI = i;
					i = (i + 1) % curve->size();
				}
			}

			prevJ = prevI;
			j = i;
			prevI = i;	// set candidate sample right after point x

			if (!isCompleted)
			{
				// bisect to find exact maximum conforming u (of midpoint between samples) - if not already |prevT,prevJ|/lfs>epsilon
				DPoint p;
				double u[2] = { (double)prevJ, (j == 0) ? (double)curve->size() : (double)j };

				while (u[1] - u[0] > EPSILON_T)
				{
					isConforming = true;
					double newU = 0.5*(u[0] + u[1]);
					p = evalCBC(&cbcs[k][4*((int)newU/sampleCount)], fmod(newU, sampleCount)/sampleCount);

					if (error > 0.0)
						isConforming = (distancePFromLine(p, prevSample, candSample) < error);

					if (isConforming)
					{
						double lfs = computeLFSForCurveAt(cbcs, *curve, medialPoints, orient[k], k, newU, ((newU - 1.0) > 0.0) ? (newU - 1.0) : 0.0, kdtree, sampleCount);
						double dist = p.distance(prevSample);
						isConforming = (dist/lfs < maxEps);
					}

					// bisect range
					if (isConforming)
						u[0] = newU;
					else
						u[1] = newU;
				}

				// fixed lfs for point u[0] farthest from prevSample
				double lfs = computeLFSForCurveAt(cbcs, *curve, medialPoints, orient[k], k, u[0], ((u[0] - 1.0) > 0.0) ? (u[0] - 1.0) : 0.0, kdtree, sampleCount);
				i = prevI;	// was conforming sample
				isConforming = true;

				// first, try to increase candSample on existing curve points
				do
				{
					if (i == intOfs)
						isCompleted = true;

					candSample = (*curve)[i];

					if (error > 0.0)
						isConforming = (distancePFromLine(p, prevSample, candSample) < error);

					if (isConforming)
					{
						double dist = p.distance(candSample);
						isConforming = (dist/lfs < maxEps);

						if (isConforming)
						{
							prevI = i;
							i++;

							if (i == (int)curve->size())
								i = 0;
						}
					}

				} while (isConforming);

				// bisect to find exact maximum conforming t for new sample
				u[0] = (double)prevI;
				u[1] = (i == 0) ? (double)curve->size() : (double)i;

				while (u[1] - u[0] > EPSILON_T)
				{
					isConforming = true;
					double newU = 0.5*(u[0] + u[1]);
					candSample = evalCBC(&cbcs[k][4*((int)newU/sampleCount)], fmod(newU, sampleCount)/sampleCount);

					if (error > 0.0)
						isConforming = (distancePFromLine(p, prevSample, candSample) < error);

					if (isConforming)
					{
						double dist = p.distance(candSample);
						isConforming = (dist/lfs < maxEps);
					}

					// bisect range
					if (isConforming)
						u[0] = newU;
					else
						u[1] = newU;
				}

				if (!isCompleted)
				{
					t = u[0];
					newSample = evalCBC(&cbcs[k][4*((int)t/sampleCount)], fmod(t, sampleCount)/sampleCount);

					// check min and max distances if given
					if (minDist != 0.0)
					{
						if (prevSample.distance(newSample) < minDist)
						{
							i = floor(t) + 1;

							if (i == (int)curve->size())
								i = 0;

							// step forward until minimum distance reached
							while (prevSample.distance(newSample) < minDist)
							{
								if (i == intOfs)
									isCompleted = true;

								newSample = evalCBC(&cbcs[k][4*(i/sampleCount)], (float)(i % sampleCount)/sampleCount);

								if (i == intOfs)
									isCompleted = true;

								prevI = i;
								i++;

								if (i == (int)curve->size())
									i = 0;
							}

							t = prevI;
						}
					}

					if (maxDist != 0.0)
					{
						if (prevSample.distance(newSample) > maxDist)
						{
							i = floor(t);

							// step backward until maximum distance reached
							while (prevSample.distance(newSample) > maxDist)
							{
								newSample = evalCBC(&cbcs[k][4*(i/sampleCount)], (float)(i % sampleCount)/sampleCount);

								prevI = i;
								i--;

								if (i == -1)
									i = (int)curve->size() - 1;
							}

							t = prevI;
						}
					}
				}
			}
		} while (!isCompleted);

		if (!isClosed)
		{
			// add new sample at end of open curve
			DPoint normal = computeNormalForCurve(1.0, &cbcs[k][4*((curve->size() - 1)/sampleCount)]);
			newSample = evalCBC(&cbcs[k][4*((curve->size() - 1)/sampleCount)], 1.0);

			if (perturb != 0.0)
			{
				// perturb new sample by up to perturb*lfs(p) along its normal in any direction
				double random = (double)rand()/RAND_MAX*2 - 1.0;
				newSample = newSample + normal*(random*perturb*lfsCurve[k][sampleCount - 1]);
			}

			CBPoint cbSample;
			cbSample.x = newSample[0];
			cbSample.y = newSample[1];
			samples[k].push_back(cbSample);
			sValues[k].push_back(t);
			noise[k].push_back(perturb*lfsCurve[k][sampleCount - 1]);
		}

	}

	annDeallocPts(ann_points);
	delete kdtree;
}

/*
 * generate curve point set from SVG curveto string
 */
void generateCurvePointsFromCBC(vector<DPoint> &cbc, vector<DPoint> &curvePoints, int sampleCount)
{
	int i, j;

	// iterate all cubic bezier curves
	for (j = 0; j < (int)cbc.size()/4; j++)
	{
		// sample the cubic bezier curve by evaluating with parameter t [0..1]
		for (i = 0; i < sampleCount; i++)
		{
			double t = (float)i/sampleCount;
			DPoint p = evalCBC(&cbc[j*4], t);
			curvePoints.push_back(p);
		}
	}
}

bool readSVGFile(string filename, vector<pair<vector<DPoint>, bool> > &curvePoints, vector<vector<DPoint> > &cbcs, list<CBEdge> &origEdges, int sampleCount)
{
	int i, pCount = 0;
	NSVGimage* nsvg = nsvgParseFromFile(filename.c_str(), "px", 96.0f);

	if (nsvg == NULL)
		return false;

	for (NSVGshape *shape = nsvg->shapes; shape != NULL; shape = shape->next)
	{
		for (NSVGpath *path = shape->paths; path != NULL; path = path->next)
		{
			curvePoints.push_back(pair<vector<DPoint>, bool>(vector<DPoint>(), path->closed));
			cbcs.push_back(vector<DPoint>());
			vector<DPoint> *cbc = &cbcs.back();

			for (i = 0; i < (path->closed ? (path->npts - 4) : (path->npts - 1)); i += 3)
			{
				float *p = &path->pts[i*2];
				cbc->push_back(DPoint(p[0], p[1]));
				cbc->push_back(DPoint(p[2], p[3]));
				cbc->push_back(DPoint(p[4], p[5]));
				cbc->push_back(DPoint(p[6], p[7]));
			}

			generateCurvePointsFromCBC(*cbc, curvePoints.back().first, sampleCount);

			for (i = 0; i < (int)curvePoints.back().first.size() - 1; i++)
			{
				CBEdge edge;
				edge.v[0] = pCount + i;
				edge.v[1] = pCount + i + 1;
				origEdges.push_back(edge);
			}

			if (path->closed)
			{
				CBEdge edge;
				edge.v[0] = pCount + curvePoints.back().first.size() - 1;
				edge.v[1] = pCount;
				origEdges.push_back(edge);
			}

			pCount += curvePoints.back().first.size();
		}
	}

	return true;
}

int generateSamplesFromSVG(string infilename, vector<vector<DPoint> > &cbcs, vector<pair<vector<DPoint>, bool> > &curvePoints,
		vector<DPoint> &medialPoints, vector<vector<bool> > &orient, vector<vector<double> > &lfs,
		vector<vector<CBPoint> > &points, vector<vector<double> > &pValues, int sampleCount,
		const float eps, float minDist, float maxDist, float &bbDiag)
{
	const float delta = 0.0;
	const float offset = 0.0;
	list<CBEdge> origEdges;
	vector<vector<float> > noise;

	readSVGFile(infilename, curvePoints, cbcs, origEdges, sampleCount);

	// determine bound box of control points
	Point min, max;
	min[0] = numeric_limits<float>::max();
	min[1] = numeric_limits<float>::max();
	max[0] = -numeric_limits<float>::max();
	max[1] = -numeric_limits<float>::max();

	for (auto cbc:cbcs)
	{
		for (auto p:cbc)
		{
			if (p[0] < min[0])
				min[0] = p[0];

			if (p[1] < min[1])
				min[1] = p[1];

			if (p[0] > max[0])
				max[0] = p[0];

			if (p[1] > max[1])
				max[1] = p[1];
		}
	}

	bbDiag = max.distance(min);

	orient.resize(curvePoints.size());
	lfs.resize(curvePoints.size());
	points.resize(curvePoints.size());
	pValues.resize(curvePoints.size());
	noise.resize(curvePoints.size());
	sampleCurvesByEps(curvePoints, cbcs, eps, 0.0, delta, points, pValues, noise, offset, medialPoints,
			orient, lfs, sampleCount, minDist*bbDiag, maxDist*bbDiag);

	return cbcs.size();
}

/*
 * write SVG file with points and edges
 */
void writeSVGFile(string filename, vector<pair<CGPoint, CGPoint> > &edgesBoundary, vector<pair<CGPoint, CGPoint> > &edgesMesh)
{
	int i;

	// constants to control SVG output can be changed here:
	const int SVG_PIXEL_DIM = 1000;	// x and y resolution on screen
	const float SVG_LINE_WIDTH1 = 1.0;	// stroke size of boundary edges
	const float SVG_LINE_WIDTH2 = 2.0;	// stroke size of boundary edges
	const float SVG_DOT_RADIUS = 6.0;	// radius of dots

	string svgHeaderPart1="<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
		"<!-- Created with CurveBenchmark -->\n"
		"<svg\n"
		"	xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n"
		"	xmlns:cc=\"http://creativecommons.org/ns#\"\n"
		"	xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n"
		"	xmlns:svg=\"http://www.w3.org/2000/svg\"\n"
		"	xmlns=\"http://www.w3.org/2000/svg\"\n"
		"	xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\"\n"
		"	xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\"\n"
		"	width=\"";
	string svgHeaderPart2="\"\n"
		"	height=\"";
	string svgHeaderPart3="\"\n"
		"	id=\"svg2\"\n"
		"	version=\"1.1\"\n"
//		"	inkscape:version=\"0.48.1 r9760\"\n"
		"	sodipodi:docname=\"";
	string svgHeaderPart4 = "\">\n"
		"<defs\n"
		"	id=\"defs4\" />\n"
		"<metadata\n"
		"	id=\"metadata7\">\n"
		"	<rdf:RDF>\n"
		"		<cc:Work\n"
		"			rdf:about=\"\">\n"
		"			<dc:format>image/svg+xml</dc:format>\n"
		"			<dc:type\n"
		"				rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />\n"
		"			<dc:title></dc:title>\n"
		"		</cc:Work>\n"
		"	</rdf:RDF>\n"
		"</metadata>\n"
		"<g\n"
		"	inkscape:label=\"Layer 1\"\n"
		"	inkscape:groupmode=\"layer\"\n"
		"	id=\"layer1\">\n";
	string svgFooter = "	</g>\n"
		"</svg>";

	// determine maximum extension in either dimension
	double min[2] = { numeric_limits<double>::max(), numeric_limits<double>::max() };
	double max[2] = { -numeric_limits<double>::max(), -numeric_limits<double>::max() };

	for (auto edge:edgesBoundary)
	{
		CGPoint p[2] = { edge.first, edge.second };

		for (i = 0; i < 2; i++)
		{
			if (p[i].x() < min[0])
				min[0] = p[i].x();

			if (p[i].x() > max[0])
				max[0] = p[i].x();

			if (p[i].y() < min[1])
				min[1] = p[i].y();

			if (p[i].y() > max[1])
				max[1] = p[i].y();
		}
	}

	for (auto edge:edgesMesh)
	{
		CGPoint p[2] = { edge.first, edge.second };

		for (i = 0; i < 2; i++)
		{
			if (p[i].x() < min[0])
				min[0] = p[i].x();

			if (p[i].x() > max[0])
				max[0] = p[i].x();

			if (p[i].y() < min[1])
				min[1] = p[i].y();

			if (p[i].y() > max[1])
				max[1] = p[i].y();
		}
	}

	double dim[2] = { max[0] - min[0], max[1] - min[1] };
	double maxDim = (dim[0] > dim[1]) ? dim[0] : dim[1];
	double factor = (SVG_PIXEL_DIM - 2*SVG_DOT_RADIUS)/maxDim;

	ofstream svgFile;
	svgFile.open(filename.data());
	svgFile << svgHeaderPart1 << SVG_PIXEL_DIM << svgHeaderPart2 << SVG_PIXEL_DIM << svgHeaderPart3 << filename.substr(0, filename.length() - 4) << svgHeaderPart4;

	// output edges as list of edge paths
	for (auto edge:edgesBoundary)
	{
		svgFile << "	<path style=\"fill:none;stroke:red;stroke-width:" << SVG_LINE_WIDTH2 << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" d=\"M";

		CGPoint p[2] = { edge.first, edge.second };

		for (i = 0; i < 2; i++)
			svgFile << " " << ((p[i].x() - min[0])*factor + SVG_DOT_RADIUS) << "," << ((p[i].y() - min[1])*factor + SVG_DOT_RADIUS);

		svgFile << " z\" id=\"boundary\" inkscape:connector-curvature=\"0\" />\n";

		for (i = 0; i < 2; i++)
			svgFile << "	<circle cx=\"" << ((p[i].x() - min[0])*factor + SVG_DOT_RADIUS) << "\" cy=\""  << ((p[i].y() - min[1])*factor + SVG_DOT_RADIUS) << "\" r=\"" << SVG_DOT_RADIUS << "\" fill=\"#0\" id=\"circle" << i << "\"/>\n";

	}

	// output edges as list of edge paths
	for (auto edge:edgesMesh)
	{
		svgFile << "	<path style=\"fill:none;stroke:lightgrey;stroke-width:" << SVG_LINE_WIDTH1 << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\" d=\"M";

		CGPoint p[2] = { edge.first, edge.second };

		for (i = 0; i < 2; i++)
			svgFile << " " << ((p[i].x() - min[0])*factor + SVG_DOT_RADIUS) << "," << ((p[i].y() - min[1])*factor + SVG_DOT_RADIUS);

		svgFile << " z\" id=\"boundary\" inkscape:connector-curvature=\"0\" />\n";
	}

	svgFile << svgFooter;
	svgFile.close();
}

void meshSampledCurve(string infilename, string outfilename, float epsilon, float minDist, float maxDist, float maxEdgeLen)
{
	const int SAMPLECOUNT = 1000;
	const int ANGLEHISTSIZE = 180;
	int i, j;
	float bbDiag;
	vector<vector<DPoint> > cbcs;
	vector<double> cbcSize;
	vector<pair<vector<DPoint>, bool> > curvePoints;
	vector<DPoint> medialPoints;
	vector<vector<bool> > orient;
	vector<vector<double> > lfsCurve, pValues;
	vector<vector<CBPoint> > points;

	PrecTimer timer;
	timer.reset();
	timer.start();
	generateSamplesFromSVG(infilename, cbcs, curvePoints, medialPoints, orient, lfsCurve, points, pValues,
			SAMPLECOUNT, epsilon, minDist, maxDist, bbDiag);
	timer.stop();
	float sampleTime = timer.getReal();

	// create conforming DT with generated polygon
	CDT cdt;

	// insert samples as vertices and edges
	Vertex_handle firstVH, prevVH, currVH;
	int pointCount = 0, edgeCount = 0;

	for (j = 0; j < (int)curvePoints.size(); j++)
	{
		for (i = 0; i < (int)points[j].size(); i++)
		{
			prevVH = currVH;
			currVH = cdt.insert(CGPoint(points[j][i].x, points[j][i].y));

			if (i == 0)
				firstVH = currVH;
			else
				if (prevVH != currVH)	// can happen at end of open curves (why?)
				{
					try
					{
						cdt.insert_constraint(prevVH, currVH);
					} catch (...)
					{
						cerr << "Exception (Intersection_of_constraints_exception?) occurred" << endl;
					}
				}
		}

		// only connect last vertex back to first vertex if curve is closed
		if (curvePoints[j].second)
			try
			{
				cdt.insert_constraint(currVH, firstVH);
			} catch (...)
			{
				cerr << "Exception (Intersection_of_constraints_exception?) occurred" << endl;
			}

		pointCount += points[j].size();
		edgeCount += points[j].size();

		if (!curvePoints[j].second)
			edgeCount--;
	}

	// add bounding box corners as samples to mesh exterior as well
	CBPoint cornerP[4];
	cornerP[0].x = numeric_limits<double>::max();
	cornerP[0].y = numeric_limits<double>::max();
	cornerP[2].x = -numeric_limits<double>::max();
	cornerP[2].y = -numeric_limits<double>::max();

	for (j = 0; j < (int)curvePoints.size(); j++)
		for (auto p:points[j])
		{
			if (p.x < cornerP[0].x)
				cornerP[0].x = p.x;

			if (p.y < cornerP[0].y)
				cornerP[0].y = p.y;

			if (p.x > cornerP[2].x)
				cornerP[2].x = p.x;

			if (p.y > cornerP[2].y)
				cornerP[2].y = p.y;
		}

#define PADDING
#ifdef PADDING
	// add 0.1 of bounding box width and height at each corner
	float width = cornerP[2].x - cornerP[0].x;
	float height = cornerP[2].y - cornerP[0].y;
	float padding = maxEdgeLen;

	if (minDist > padding)
		padding = minDist;

	if (maxDist > padding)
		padding = maxDist;

	cornerP[0].x -= padding*width;
	cornerP[2].x += padding*width;
	cornerP[0].y -= padding*height;
	cornerP[2].y += padding*height;

	cornerP[1].x = cornerP[2].x;
	cornerP[1].y = cornerP[0].y;
	cornerP[3].x = cornerP[0].x;
	cornerP[3].y = cornerP[2].y;

	Vertex_handle cornerVH[4];

	for (i = 0; i < 4; i++)
		cornerVH[i] = cdt.insert(CGPoint(cornerP[i].x, cornerP[i].y));

	for (i = 0; i < 4; i++)
		cdt.insert_constraint(cornerVH[i], cornerVH[(i + 1) % 4]);
#endif

	cout << "number of samples: " << pointCount << ", vertices on boundary: " << cdt.number_of_vertices() << endl;

	Mesher mesher(cdt);

	timer.reset();
	timer.start();
	mesher.set_criteria(Criteria(0.125, maxEdgeLen*bbDiag));	// arcsin(1/(2*sqrt(1/(4*b))))*180/pi= (1/8 -> 20.7, 1/2 -> 45)
	mesher.refine_mesh();
	lloyd_optimize_mesh_2(cdt, CGAL::parameters::max_iteration_number = 10);
	timer.stop();
	float meshingTime = timer.getReal();

	cout << "number of vertices in mesh: " << cdt.number_of_vertices() << ", faces: " << cdt.number_of_faces() << endl;

	// compute angle histogram
	int angles[ANGLEHISTSIZE];

	for (i = 0; i < ANGLEHISTSIZE; i++)
		angles[i] = 0;

	float minAngle = PI, maxAngle = 0.0;
	int triangleCount = 0;

	for (Triangulation_2::Finite_faces_iterator iter = cdt.finite_faces_begin(); iter != cdt.finite_faces_end(); iter++)
		if (iter->is_in_domain())
		{
			// compute triangle's angles
			CGPoint p[3];
			Point v[3];

			// triangle's points
			for (i = 0; i < 3; i++)
				p[i] = iter->vertex(i)->point();

			// triangle's edge vectors
			for (i = 0; i < 3; i++)
			{
				v[i] = Point(p[(i + 1) % 3].x() - p[i].x(), p[(i + 1) % 3].y() - p[i].y());
				v[i].normalize();
			}

			// compute angles
			for (i = 0; i < 3; i++)
			{
				Point prevV(-v[(i + 2) % 3].x(), -v[(i + 2) % 3].y());
				float angle = acos(v[i].dot(prevV));

				if (angle < minAngle)
					minAngle = angle;

				if (angle > maxAngle)
					maxAngle = angle;

				angles[(int)(angle/PI*ANGLEHISTSIZE)]++;
			}

			triangleCount++;
		}

	minAngle *=180.0/PI;
	maxAngle *=180.0/PI;

	cout << "csv: name,vertices,triangles,sampletime,meshingtime" << endl;
	cout << "csv," << outfilename << "," << pointCount << "," << triangleCount << "," << sampleTime << "," << meshingTime << endl;

	ofstream outfile;

	string stlfilename = outfilename + ".stl";
	string jigfilename = outfilename + "_jigsaw.msh";
	string plyfilename = outfilename + ".ply";

	outfile.open(stlfilename);
    outfile << "solid name\n";

	for (Triangulation_2::Finite_faces_iterator iter = cdt.finite_faces_begin(); iter != cdt.finite_faces_end(); iter++)
		if (iter->is_in_domain())
		{
			outfile << "facet normal 0.0 0.0 0.0\nouter loop\n";

			for (i = 0; i < 3; i++)
				outfile << "vertex " << iter->vertex(i)->point() << " 0.0\n";

			outfile << "endloop\nendfacet\n";
		}

    outfile << "endsolid name\n";
	outfile.close();

	// write polylines between samples to PLY file
	outfile.open(plyfilename);
	outfile << "ply\nformat ascii 1.0\nelement vertex " << pointCount << "\nproperty float x\nproperty float y\nproperty float z\n" <<
			"element face 0\nproperty list uchar int vertex_indices\n" <<
			"element edge " << edgeCount << "\nproperty int vertex1\nproperty int vertex2\nend_header" << endl;

	for (j = 0; j < (int)curvePoints.size(); j++)
		for (auto p:points[j])
			outfile << p.x << " " << p.y << " 0" << endl;

	int ofs = 0;

	for (j = 0; j < (int)curvePoints.size(); j++)
	{
		// only connect last vertex back to first vertex if curve is closed
		for (i = 0; i < (int)points[j].size() - (curvePoints[j].second ? 0 : 1); i++)
			outfile << ofs + i << " " << ofs + ((i + 1) % points[j].size()) << endl;

		ofs += points[j].size();
	}

	outfile.close();

	// determine edges between samples and remaining edges in mesh
	set<CGPoint> pointSet;

	for (j = 0; j < (int)curvePoints.size(); j++)
	{
		for (i = 0; i < (int)points[j].size(); i++)
		{
			CGPoint p(points[j][i].x, points[j][i].y);
			pointSet.insert(p);
		}
	}

	vector<pair<CGPoint, CGPoint> > edgesBoundary, edgesMesh, edgesEmpty;

	for (Triangulation_2::Finite_edges_iterator iter = cdt.finite_edges_begin(); iter != cdt.finite_edges_end(); iter++)
	{
		CGPoint p0 = iter->first->vertex(cdt.cw(iter->second))->point();
		CGPoint p1 = iter->first->vertex(cdt.ccw(iter->second))->point();

		if (cdt.is_constrained(*iter))
			edgesBoundary.push_back(pair<CGPoint, CGPoint>(p0, p1));
		else
			edgesMesh.push_back(pair<CGPoint, CGPoint>(p0, p1));
	}

	// write polylines between samples colored red to SVG file
	string svg1filename = outfilename + "-1.svg";
	writeSVGFile(svg1filename, edgesBoundary, edgesEmpty);

	// write mesh with polylines between samples colored red to SVG file
	string svg2filename = outfilename + "-2.svg";
	writeSVGFile(svg2filename, edgesBoundary, edgesMesh);
}

int main(int argc, char **argv)
{
	float minDist = 0.0, maxDist = 0.0, maxEdgeLen = 0.0;

	if (argc < 4)
	{
		cerr << "sample2d input.svg outfile-basename epsilon [mindist] [maxdist] [maxedgelen]\nall optional arguments are in function of the bounding box diagonal of the control points" << endl;
		exit(1);
	}

	if (argc > 4)
		minDist = atof(argv[4]);

	if ((minDist < 0.0) || (minDist >= 1.0))
	{
		cerr << "minDist must be in range [0,1[ (of bounding box diagonal, 0 == no minDist)" << endl;
		exit(1);
	}

	if (argc > 5)
		maxDist = atof(argv[5]);

	if ((maxDist < 0.0) || (maxDist >= 1.0))
	{
		cerr << "maxDist must be in range [0,1[ (of bounding box diagonal, 0 == no maxDist)" << endl;
		exit(1);
	}

	if (argc > 6)
		maxEdgeLen = atof(argv[6]);

	if ((maxEdgeLen < 0.0) || (maxEdgeLen >= 1.0))
	{
		cerr << "maxEdgeLen must be in range [0,1[ (of bounding box diagonal, 0 == no maxEdgeLen)" << endl;
		exit(1);
	}

	if ((minDist > 0) && (maxDist > 0))
		assert(minDist <= maxDist);

	meshSampledCurve(string(argv[1]), string(argv[2]), atof(argv[3]), minDist, maxDist, maxEdgeLen);
}
