#pragma once 
#include <geos.h>
#include <float.h>
#include <math.h>
using namespace std;

struct WayCoverLength
{
	double len_ln1_cover;
	double len_ln2_cover;
	WayCoverLength()
	 :len_ln1_cover(0.)
	 ,len_ln2_cover(0.){}
	~WayCoverLength(){}
};

struct WayMappingPt
{
	geos::geom::Coordinate startp;
	geos::geom::Coordinate endp;
	int startp_pos;
	int endp_pos;
	geos::geom::Coordinate startp_mappingp;
	geos::geom::Coordinate endp_mappingp;
	int startp_mappingp_pos;
	int endp_mappingp_pos;
	WayMappingPt()
	 :startp(geos::geom::Coordinate())
	 ,endp(geos::geom::Coordinate())
	 ,startp_pos(-1)
	 ,endp_pos(-1)
	 ,startp_mappingp(geos::geom::Coordinate())
	 ,endp_mappingp(geos::geom::Coordinate())
	 ,startp_mappingp_pos(-1)
	 ,endp_mappingp_pos(-1){}
	~WayMappingPt(){}
};

struct WayMappingPos
{
	WayMappingPt ln1;
	WayMappingPt ln2;
	bool mapping_ln1;
	bool mapping_ln2;
	WayMappingPos()
	 :ln1(WayMappingPt())
	 ,ln2(WayMappingPt())
	 ,mapping_ln1(false)
	 ,mapping_ln2(false){}
	~WayMappingPos(){}
};

struct WayMappingDis
{
	double dis_min;
	double dis_max;
	double dis_avg;
	WayMappingDis()
	 : dis_min(DBL_MAX)
	 , dis_max(0.)
	 , dis_avg(DBL_MAX){}

	~WayMappingDis(){}

};

struct WayHaustorffDis
{
	double dis_hausdorff_ln1;
	double dis_hausdorff_ln2;
	WayHaustorffDis()
	 : dis_hausdorff_ln1(DBL_MAX)
	 , dis_hausdorff_ln2(DBL_MAX){}

	~WayHaustorffDis(){}
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

	double dis_hausdorff_ln1;
	double dis_hausdorff_ln2;
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
	, dis_hausdorff_ln1(DBL_MAX)
	, dis_hausdorff_ln2(-1.)
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

	// 求sg2 p0点作sg2的垂线，与sg1的交点 sg1与sg2的夹角小于30度
	static void calMappingPoint(geos::geom::Coordinate& pt, const geos::geom::LineSegment& sg1, const geos::geom::LineSegment& sg2, bool first);

	//	ln1映射到ln2的起止坐标点和对应的位置
	static bool calMappingPointPos(WayMappingPt& wayMappingPt, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

	// ln1与ln2, 相互到第一个匹配点和最后一个匹配点，以及对应到映射点
	static void calCoverOfLine(WayMappingPos& wayMappingPos, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

	// ln1与ln2之间点豪斯多夫距离 包括ln1到ln2和ln2到ln1的两个豪斯多夫距离的最大和最小值
	static bool calHausdorffDis(WayHaustorffDis& wayHausdorffDis, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2);

	// ln1与ln2匹配点间的距离 包括最大／最小／平均
	static bool calMappingDis(WayMappingDis& wayMappingDis, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

	// 计算ln1于ln2分别的匹配长度
	static void calMappingLength(WayCoverLength& wayCoverLength, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

	// ln1和ln2之间的豪斯多夫距离和匹配距离
	static bool calMappingFeature(WayFeature& wayFeature, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);
};
