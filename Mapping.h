#pragma once 
#include <geos.h>
#include <float.h>
#include <math.h>
using namespace std;

struct CoverInfo
{
	geos::geom::Coordinate startp_ln1;
	geos::geom::Coordinate endp_ln1;
	int startp_ln1_pos;
	int endp_ln1_pos;
	geos::geom::Coordinate startp_ln1_mappingp;
	geos::geom::Coordinate endp_ln1_mappingp;
	int startp_ln1_mappingp_pos;
	int endp_ln1_mappingp_pos;

	geos::geom::Coordinate startp_ln2;
	geos::geom::Coordinate endp_ln2;
	int startp_ln2_pos;
	int endp_ln2_pos;
	geos::geom::Coordinate startp_ln2_mappingp;
	geos::geom::Coordinate endp_ln2_mappingp;
	int startp_ln2_mappingp_pos;
	int endp_ln2_mappingp_pos;

	bool covered = true;
};

//特征计算（平均匹配距离，最小距离，最大距离，目标覆盖率，参考覆盖率，匹配角度差，进入角度差，离开角度差，首尾角度差，内角和差，豪斯多夫距离，相交点数，入度，出度）
struct WayFeature
{
	double cover_ln1_ln2;
	double cover_ln2_ln1;

	double angle_diff_mapping;
	double angle_diff_start;
	double angle_diff_end;
	double angle_diff_skeleton;
	double angle_diff_all;

	double dis_hausdorff_min;
	double dis_hausdorff_max;
	double dis_hausdorff_avg;
	double dis_average;
	double dis_max;
	double dis_min;

	int insection_point_num;

	int in_degree;
	int out_degree;

	double length_mapping_ln1;
	double length_mapping_ln2;

	double length_ln1;
	double length_ln2;

	WayFeature()
	: cover_ln1_ln2(0.)
	, cover_ln2_ln1(0.)
	, angle_diff_mapping(M_PI)
	, angle_diff_start(M_PI)
	, angle_diff_end(M_PI)
	, angle_diff_skeleton(M_PI)
	, angle_diff_all(M_PI)
	, dis_hausdorff_min(DBL_MAX)
	, dis_hausdorff_max(-1.)
	, dis_hausdorff_avg(DBL_MAX)
	, dis_average(DBL_MAX)
	, dis_max(-1.)
	, dis_min(DBL_MAX)
	, insection_point_num(0)
	, in_degree(0)
	, out_degree(0)
	, length_mapping_ln1(0.)
	, length_mapping_ln2(0.)
	, length_ln1(0)
	, length_ln2(0){}
	~WayFeature() {}
};

class Mapping{
	
public:
	Mapping();
	virtual ~Mapping();
	
//	static double pt2Segment(geos::geom::Point& p, geos::geom::LineString& sg);
	static geos::geom::Coordinate calPt_pt2Segment(geos::geom::Point& p, geos::geom::Point& p1, geos::geom::Point& p2);
	static geos::geom::Coordinate calPt_pt2Segment(geos::geom::Point* p, geos::geom::Point* p1, geos::geom::Point* p2);
//	static double pt2Line(geos::geom::Point p, geos::geom::LineString& ls);
	static geos::geom::Coordinate calPt_Pt2Line(geos::geom::Point& p, geos::geom::LineString& ls);
	static geos::geom::Coordinate calPt_Pt2Line(geos::geom::Point* p, geos::geom::LineString* ls);

	static CoverInfo cal_CoverOfLine(geos::geom::LineString& ln1, geos::geom::LineString& ln2, double dist_Threshold);

	static geos::geom::Coordinate calMappingPoint(geos::geom::LineSegment& sg1, geos::geom::LineSegment& sg2);

	static WayFeature calMappingFeature(geos::geom::LineString& ln1, geos::geom::LineString& ln2, double dis_Threshold);
};
