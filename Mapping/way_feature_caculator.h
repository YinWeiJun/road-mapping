#pragma once 
#include "../common/geos_tool.h"
#include <geos.h>
#include <float.h>
#include <math.h>
using namespace std;

#include "way_matching_operator.h"

//#include <geos/geom/LineString.h>

struct WayCoverLength
{
  double len_ln1_cover;
  double len_ln2_cover;
  double angle_diff_cover;
  WayCoverLength()
   : len_ln1_cover(0.)
   , len_ln2_cover(0.)
   , angle_diff_cover(DBL_MAX){}
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
	double dis_mapping_sta;
	double dis_mapping_end;
	WayMappingDis()
	 : dis_min(DBL_MAX)
	 , dis_max(0.)
	 , dis_avg(DBL_MAX)
	 , dis_mapping_sta(DBL_MAX)
	 , dis_mapping_end(DBL_MAX){}

	~WayMappingDis(){}

};

struct WayHausdorffDis
{
	double dis_hausdorff_ln1;
	double dis_hausdorff_ln2;
	WayHausdorffDis()
	 : dis_hausdorff_ln1(DBL_MAX)
	 , dis_hausdorff_ln2(DBL_MAX){}

	~WayHausdorffDis(){}
};

struct WayHausdorffDir
{
	double dir_hausdorff_ln1;
	double dir_hausdorff_ln2;
	WayHausdorffDir()
		: dir_hausdorff_ln1(DBL_MAX)
		, dir_hausdorff_ln2(DBL_MAX) {}

	~WayHausdorffDir(){}
};

class WayFeatureCalculator
{
public:
	WayFeatureCalculator();
	virtual ~WayFeatureCalculator();

	// 求sg2 p0点作sg2的垂线，与sg1的交点 sg1与sg2的夹角小于30度
	static void calMappingPoint(geos::geom::Coordinate& pt, const geos::geom::LineSegment& sg1, const geos::geom::LineSegment& sg2, bool first);

	//	ln1映射到ln2的起止坐标点和对应的位置
	static bool calMappingPointPos(WayMappingPt& wayMappingPt, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

	// 计算ln1于ln2分别的匹配长度
	static bool calMappingLength(WayCoverLength& wayCoverLength, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

	// ln1与ln2, 相互到第一个匹配点和最后一个匹配点，以及对应到映射点
	static bool calCoverOfLine(WayMappingPos& wayMappingPos, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

	// ln1与ln2之间点豪斯多夫距离 包括ln1到ln2和ln2到ln1的两个豪斯多夫距离的最大和最小值
	static bool calHausdorffDis(WayHausdorffDis& wayHausdorffDis, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2);

	// ln1与ln2匹配点间的距离 包括最大／最小／平均
	static bool calMappingDis(WayMappingDis& wayMappingDis, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

	// 计算ln1与ln2的豪斯多夫角度
	static bool calHausdorffDir(WayHausdorffDir& wayHausdorffDir, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2);

	// ln1和ln2之间的豪斯多夫距离和匹配距离
	static bool calMappingFeature(WayFeature& wayFeature, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold);

};
