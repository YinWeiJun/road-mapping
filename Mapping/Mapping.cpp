/*
 * Mapping.cpp
 *
 *  Created on: 2018年5月4日
 *      Author: didi
 */
#include <geos.h>
#include "Mapping.h"
#include <math.h>
using namespace std;
const double ANGLE_THRESHOLD = M_PI/6;

Mapping::Mapping()
{
}

Mapping::~Mapping()
{
}

// 求sg2 p0点作sg2的垂线，与sg1的交点 sg1与sg2的夹角小于30度
void Mapping::calMappingPoint(geos::geom::Coordinate& insp, const geos::geom::LineSegment& sg1, const geos::geom::LineSegment& sg2, bool first)
{
	double k1 = 0;
	double b1 = 0;
	double k2 = 0;
	double b2 = 0;

	geos::geom::Coordinate pt;
	if(first)
		pt = sg2.p0;
	else
		pt = sg2.p1;
	if(sg2.p0.y == sg2.p1.y)
	{
		insp.x = pt.x;
		if(sg1.p0.y == sg1.p1.y)
		{
			insp.y = sg1.p0.y;
		}
		else
		{
			k1 = (sg1.p1.y - sg1.p0.y) / (sg1.p1.x - sg1.p0.x);
			b1 = sg1.p0.y - k1 * sg1.p0.x;
			insp.y = k1 * insp.x + b1;
		}
	}
	else if(sg2.p0.x == sg2.p1.x)
	{
		insp.y = pt.y;
		if(sg1.p0.x == sg1.p1.x)
		{
			insp.x = sg1.p0.x;
		}
		else
		{
			k1 = (sg1.p1.y - sg1.p0.y) / (sg1.p1.x - sg1.p0.x);
			b1 = sg1.p0.y - k1 * sg1.p0.x;
			insp.x = (insp.y - b1) / k1;
		}
	}
	else
	{
		k2 = (sg2.p1.y - sg2.p0.y) / (sg2.p1.x - sg2.p0.x);
		k2 = -1.0 / k2;
		b2 = sg2.p0.y - k2 * sg2.p0.x;
		if(sg1.p0.x == sg1.p1.x)
		{
			insp.x = pt.x;
			insp.y = k2 * insp.x + b2;
		}
		else
		{
			k1 = (sg1.p1.y - sg1.p0.y) / (sg1.p1.x - sg1.p0.x);
			b1 = sg1.p0.y - k1 * sg1.p0.x;
			insp.x = (b2 - b1) / (k1 - k2);
			insp.y = k1 * insp.x + b1;
		}
	}
}

// ln1映射到ln2的起止坐标点和对应的位置
bool Mapping::calMappingPointPos(WayMappingPt& wayMappingPt, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	double angle_Threshold = ANGLE_THRESHOLD;
	bool finded = false;
	for(int i = 1; i<ln1.getNumPoints() && !finded; i++)
	{
		geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(i-1), ln1.getCoordinateN(i));
		geos::geom::Coordinate tmp = ln1.getCoordinateN(i-1);
		double angle_sg_ln1 = sg_ln1.angle();
		for(int j = 1; j<ln2.getNumPoints(); j++)
		{
			geos::geom::Coordinate prjc = geos::geom::Coordinate();
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
			double angle_sg_ln2 = sg_ln2.angle();
			sg_ln2.project(tmp, prjc);
			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
			double d_angle2 = 0;
			if (i>1)
			{
				geos::geom::LineSegment sg_ln1_bef = geos::geom::LineSegment(ln1.getCoordinateN(i-2), ln1.getCoordinateN(i-1));
				d_angle2 = fabs(sg_ln1_bef.angle() - angle_sg_ln2);
			}
			if (sg_ln2.distance(prjc) < 1e-7 && tmp.distance(prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < angle_Threshold))
			{
				wayMappingPt.startp = tmp;
				wayMappingPt.startp_pos = i-1;
				wayMappingPt.startp_mappingp = prjc;
				if (prjc.equals(sg_ln2.p0))
					wayMappingPt.startp_mappingp_pos = j - 1;
				else
					wayMappingPt.startp_mappingp_pos = j;
				if(1 < i and !prjc.equals(sg_ln2.p0))
				{
					geos::geom::LineSegment sg_ln1_tmp = geos::geom::LineSegment(ln1.getCoordinateN(i-2), ln1.getCoordinateN(i-1));
					geos::geom::Coordinate startp;
					Mapping::calMappingPoint(startp, sg_ln1_tmp, sg_ln2, true);
					if(sg_ln1_tmp.distance(startp) < 1e-7 && sg_ln2.p0.distance(startp) < dis_Threshold)
					{
						wayMappingPt.startp = startp;
						wayMappingPt.startp_pos = i-1;
						wayMappingPt.startp_mappingp = sg_ln2.p0;
						wayMappingPt.startp_mappingp_pos = j-1;
					}
				}
				finded = true;
				break;
			}
		}
		if(i == ln1.getNumPoints()-1)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(0), ln2.getCoordinateN(1));
			geos::geom::Coordinate startp;
			Mapping::calMappingPoint(startp, sg_ln1, sg_ln2, true);
			if(sg_ln1.distance(startp) < 1e-7 && sg_ln2.p0.distance(startp) < dis_Threshold)
			{
				wayMappingPt.startp = startp;
				wayMappingPt.startp_pos = i;
				wayMappingPt.startp_mappingp = sg_ln2.p0;
				wayMappingPt.startp_mappingp_pos = 0;
				finded = true;
				break;
			}
		}
	}
	bool finded2 = finded;
	for(int i = ln1.getNumPoints() - 1; i > 0 && finded2; i--)
	{
		geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(i), ln1.getCoordinateN(i-1));
		geos::geom::Coordinate tmp = ln1.getCoordinateN(i);
		double angle_sg_ln1 = sg_ln1.angle();
		for(int j = ln2.getNumPoints() - 1; j > 0; j--)
		{
			geos::geom::Coordinate prjc = geos::geom::Coordinate();
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j), ln2.getCoordinateN(j-1));
			double angle_sg_ln2 = sg_ln2.angle();
			sg_ln2.project(tmp, prjc);
			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
			double d_angle2 = 0;
			if (i < ln2.getNumPoints() - 1)
			{
				geos::geom::LineSegment sg_ln1_bef = geos::geom::LineSegment(ln1.getCoordinateN(i+1), ln1.getCoordinateN(i));
				d_angle2 = fabs(sg_ln1_bef.angle() - angle_sg_ln2);
			}
			if (sg_ln2.distance(prjc) < 1e-7 && tmp.distance(prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < angle_Threshold))
			{
				wayMappingPt.endp = tmp;
				wayMappingPt.endp_pos = i;
				wayMappingPt.endp_mappingp = prjc;
				if (prjc.equals(sg_ln2.p0))
					wayMappingPt.endp_mappingp_pos = j;
				else
					wayMappingPt.endp_mappingp_pos = j - 1;
				if(i < ln1.getNumPoints() - 1 and !prjc.equals(sg_ln2.p0))
				{
					geos::geom::LineSegment sg_ln1_tmp = geos::geom::LineSegment(ln1.getCoordinateN(i+1), ln1.getCoordinateN(i));
					geos::geom::Coordinate startp;
					Mapping::calMappingPoint(startp, sg_ln1_tmp, sg_ln2, true);
					if(sg_ln1_tmp.distance(startp) < 1e-7 && sg_ln2.p0.distance(startp) < dis_Threshold)
					{
						wayMappingPt.endp = startp;
						wayMappingPt.endp_pos = i;
						wayMappingPt.endp_mappingp = sg_ln2.p0;
						wayMappingPt.endp_mappingp_pos = j;
					}
				}
				finded2 = false;
				break;
			}
		}
		if(1 == i)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(ln2.getNumPoints()-1), ln2.getCoordinateN(ln2.getNumPoints()-1));
			geos::geom::Coordinate endp;
			Mapping::calMappingPoint(endp, sg_ln1, sg_ln2, dis_Threshold);
			if(sg_ln1.distance(endp) < 1e-7 && sg_ln2.p0.distance(endp) < dis_Threshold)
			{
				wayMappingPt.endp = endp;
				wayMappingPt.endp_pos = 0;
				wayMappingPt.endp_mappingp = sg_ln2.p0;
				wayMappingPt.endp_mappingp_pos = ln2.getNumPoints() -1;
				finded2 = false;
				break;
			}
		}
	}

	return finded;
}

// ln1与ln2, 相互到第一个匹配点和最后一个匹配点，以及对应到映射点
void Mapping::calCoverOfLine(WayMappingPos& wayMappingPos, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	WayMappingPt wayMappingPt_ln1;
	bool mapping_ln1 = Mapping::calMappingPointPos(wayMappingPt_ln1, ln1, ln2, dis_Threshold);
	if(mapping_ln1)
	{
		wayMappingPos.ln1 = wayMappingPt_ln1;
	}
	wayMappingPos.mapping_ln1 = mapping_ln1;

	WayMappingPt wayMappingPt_ln2;
	bool mapping_ln2 = Mapping::calMappingPointPos(wayMappingPt_ln2, ln2, ln1, dis_Threshold);
	if(mapping_ln2)
	{
		wayMappingPos.ln2 = wayMappingPt_ln2;
	}
	wayMappingPos.mapping_ln2 = mapping_ln2;

}

// ln1与ln2之间点豪斯多夫距离 包括ln1到ln2和ln2到ln1的两个豪斯多夫距离的最大和最小值
bool Mapping::calHausdorffDis(WayHaustorffDis& wayHausdorff, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2)
{
	double dis_hausdorff_ln1 = 0;
	int j = 1;
	for(int i=0; i<ln1.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		double d = ln2.getCoordinateN(j-1).distance(pt_ln1);
		for( ; j<ln2.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
			double dis = sg_ln2.distance(pt_ln1);
			if(dis < d)
				d = dis;
			else
				break;
		}
		dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, d);
	}

	double dis_hausdorff_ln2 = 0;
	j = 1;
	for(int i=0; i<ln2.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln2 = ln2.getCoordinateN(i);
		double d = ln1.getCoordinateN(j-1).distance(pt_ln2);
		for( ; j<ln1.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
			double dis = sg_ln1.distance(pt_ln2);
			if(dis < d)
				d = dis;
			else
				break;
		}
		dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, d);
	}
	wayHausdorff.dis_hausdorff_ln1 = dis_hausdorff_ln1;
	wayHausdorff.dis_hausdorff_ln2 = dis_hausdorff_ln2;
	return true;
}

// ln1与ln2匹配点间的距离 包括最大／最小／平均
bool Mapping::calMappingDis(WayMappingDis& wayMappingDis, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	WayMappingPos wayMappingPos;
	Mapping::calCoverOfLine(wayMappingPos, ln1, ln2, dis_Threshold);
	int num = 0;
	double dis_avg = 0;
	double dis_min = DBL_MAX;
	double dis_max = 0;
	int j = wayMappingPos.ln1.startp_mappingp_pos;
	if(0 == j)
		j += 1;
	if(!(wayMappingPos.ln1.startp.equals(ln1.getCoordinateN(wayMappingPos.ln1.startp_pos))))
	{
		double d = wayMappingPos.ln1.startp.distance(wayMappingPos.ln1.startp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}
	for(int i=wayMappingPos.ln1.startp_pos; i<=wayMappingPos.ln1.endp_pos && i<ln1.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		for( ; j<=wayMappingPos.ln1.endp_mappingp_pos && j < ln2.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
			geos::geom::Coordinate prjp = geos::geom::Coordinate();
			sg_ln2.project(pt_ln1, prjp);
			if(sg_ln2.distance(prjp) < 1e-7)
			{
				num += 1;
				double d = pt_ln1.distance(prjp);
				dis_avg += d;
				dis_min = std::min(dis_min, d);
				dis_max = std::max(dis_max, d);
				break;
			}
		}
	}
	if(!(wayMappingPos.ln1.endp.equals(ln1.getCoordinateN(wayMappingPos.ln1.endp_pos))))
	{
		double d = wayMappingPos.ln1.endp.distance(wayMappingPos.ln1.endp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}

	j = wayMappingPos.ln2.startp_mappingp_pos;
	if(0 == j)
		j += 1;
	if(!(wayMappingPos.ln2.startp.equals(ln2.getCoordinateN(wayMappingPos.ln2.startp_pos))))
	{
		double d = wayMappingPos.ln2.startp.distance(wayMappingPos.ln2.startp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}
	for(int i=wayMappingPos.ln2.startp_pos; i<=wayMappingPos.ln2.endp_pos && i<ln2.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln2 = ln2.getCoordinateN(i);
		for( ; j<=wayMappingPos.ln2.endp_mappingp_pos && j<ln1.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
			geos::geom::Coordinate prjp = geos::geom::Coordinate();
			sg_ln1.project(pt_ln2, prjp);
			if(sg_ln1.distance(prjp) < 1e-7)
			{
				num += 1;
				double d = pt_ln2.distance(prjp);
				dis_avg += d;
				dis_min = std::min(dis_min, d);
				dis_max = std::max(dis_max, d);
				break;
			}
		}
	}
	if(!(wayMappingPos.ln2.endp.equals(ln2.getCoordinateN(wayMappingPos.ln2.endp_pos))))
	{
		double d = wayMappingPos.ln2.endp.distance(wayMappingPos.ln2.endp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}

	wayMappingDis.dis_avg = dis_avg / num;
	wayMappingDis.dis_min = dis_min;
	wayMappingDis.dis_max = dis_max;
	return true;
}

// 计算ln1于ln2分别的匹配长度
void Mapping::calMappingLength(WayCoverLength& wayCoverLength, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	WayMappingPos wayMappingPos;
	Mapping::calCoverOfLine(wayMappingPos, ln1, ln2, dis_Threshold);
	double len_ln1_cover = 0.;
	int left = wayMappingPos.ln1.startp_pos;
	int right = wayMappingPos.ln1.endp_pos;
	if(!wayMappingPos.ln1.startp.equals(ln1.getCoordinateN(wayMappingPos.ln1.startp_pos)))
	{
		len_ln1_cover += (wayMappingPos.ln1.startp).distance(ln1.getCoordinateN(wayMappingPos.ln1.startp_pos));
	}
	while(left < right)
	{
		len_ln1_cover += ln1.getCoordinateN(left).distance(ln1.getCoordinateN(left+1));
		left += 1;
	}
	if(!wayMappingPos.ln1.endp.equals(ln1.getCoordinateN(wayMappingPos.ln1.endp_pos)))
	{
		len_ln1_cover += wayMappingPos.ln1.endp.distance(ln1.getCoordinateN(wayMappingPos.ln1.endp_pos));
	}

	double len_ln2_cover = 0.;
	left = wayMappingPos.ln2.startp_pos;
	right = wayMappingPos.ln2.endp_pos;
	if(!wayMappingPos.ln2.startp.equals(ln2.getCoordinateN(wayMappingPos.ln2.startp_pos)))
	{
		len_ln2_cover += (wayMappingPos.ln2.startp).distance(ln2.getCoordinateN(wayMappingPos.ln2.startp_pos));
	}
	while(left < right)
	{
		len_ln2_cover += (ln2.getCoordinateN(left)).distance(ln2.getCoordinateN(left+1));
		left += 1;
	}
	if(!wayMappingPos.ln2.endp.equals(ln2.getCoordinateN(wayMappingPos.ln2.endp_pos)))
	{
		len_ln2_cover += wayMappingPos.ln2.endp.distance(ln2.getCoordinateN(wayMappingPos.ln2.endp_pos));
	}

	wayCoverLength.len_ln1_cover = len_ln1_cover;
	wayCoverLength.len_ln2_cover = len_ln2_cover;

}

// ln1和ln2之间的豪斯多夫距离和匹配距离
bool Mapping::calMappingFeature(WayFeature& wayFeature, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	WayMappingPos wayMappingPos;
	Mapping::calCoverOfLine(wayMappingPos, ln1, ln2, dis_Threshold);
	int num = 0;
	double dis_avg = 0;
	double dis_min = DBL_MAX;
	double dis_max = 0;
	double dis_hausdorff_ln1 = 0;
	double dis_hausdorff_ln2 = 0;
	// 计算ln1起点到匹配起点间到Hausdorff距离
	for(int i=0; i<wayMappingPos.ln1.startp_pos; i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		if(wayMappingPos.ln1.startp_mappingp_pos == 0)
		{
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, pt_ln1.distance(ln2.getCoordinateN(0)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=1; j<=wayMappingPos.ln1.startp_mappingp_pos && j<ln1.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
				d = std::min(d, sg_ln2.distance(pt_ln1));
			}
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, d);
		}
	}
	// 计算ln1第一个匹配点到最后一个匹配点间的Hausdorff距离和匹配距离
	int j = wayMappingPos.ln1.startp_mappingp_pos;
	if(0 == j)
		j += 1;
	if(!(wayMappingPos.ln1.startp.equals(ln1.getCoordinateN(wayMappingPos.ln1.startp_pos))))
	{
		double d = wayMappingPos.ln1.startp.distance(wayMappingPos.ln1.startp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);

	}
	for(int i=wayMappingPos.ln1.startp_pos; i<=wayMappingPos.ln1.endp_pos && i<ln1.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		for( ; j<=wayMappingPos.ln1.endp_mappingp_pos && j<ln2.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
			geos::geom::Coordinate prjp = geos::geom::Coordinate();
			sg_ln2.project(pt_ln1, prjp);
			if(sg_ln2.distance(prjp) < 1e-7)
			{
				num += 1;
				double d = pt_ln1.distance(prjp);
				dis_avg += d;
				dis_min = std::min(dis_min, d);
				dis_max = std::max(dis_max, d);
				dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, d);
				break;
			}
		}
	}
	if(!(wayMappingPos.ln1.endp.equals(ln1.getCoordinateN(wayMappingPos.ln1.endp_pos))))
	{
		double d = wayMappingPos.ln1.endp.distance(wayMappingPos.ln1.endp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}
	// 计算最后一个匹配点到最后一个点间到Hausdorff距离
	for(int i=wayMappingPos.ln1.endp_pos+1; i<ln1.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		if(wayMappingPos.ln1.endp_mappingp_pos == ln2.getNumPoints()-1)
		{
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, pt_ln1.distance(ln2.getCoordinateN(ln2.getNumPoints()-1)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=wayMappingPos.ln1.endp_mappingp_pos+1; j<ln2.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
				d = std::min(d, sg_ln2.distance(pt_ln1));
			}
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, d);
		}
	}

	// 同上，ln2到ln1
	j = wayMappingPos.ln2.startp_mappingp_pos;
	if(0 == j)
		j += 1;
	for(int i=0; i<wayMappingPos.ln2.startp_pos; i++)
	{
		geos::geom::Coordinate pt_ln2 = ln1.getCoordinateN(i);
		if(wayMappingPos.ln2.startp_mappingp_pos == 0)
		{
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, pt_ln2.distance(ln1.getCoordinateN(0)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=1; j<=wayMappingPos.ln2.startp_mappingp_pos && j<ln2.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
				d = std::min(d, sg_ln1.distance(pt_ln2));
			}
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, d);
		}
	}

	if(!(wayMappingPos.ln2.startp.equals(ln2.getCoordinateN(wayMappingPos.ln2.startp_pos))))
	{
		double d = wayMappingPos.ln2.startp.distance(wayMappingPos.ln2.startp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
		std::cout << "--" + std::to_string(d) << std::endl;
	}
	for(int i=wayMappingPos.ln2.startp_pos; i<=wayMappingPos.ln2.endp_pos && i<ln2.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln2 = ln2.getCoordinateN(i);
		for( ; j<=wayMappingPos.ln2.endp_mappingp_pos && j<ln1.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
			geos::geom::Coordinate prjp = geos::geom::Coordinate();
			sg_ln1.project(pt_ln2, prjp);
			if(sg_ln1.distance(prjp) < 1e-7)
			{
				num += 1;
				double d = pt_ln2.distance(prjp);
				dis_avg += d;
				dis_min = std::min(dis_min, d);
				dis_max = std::max(dis_max, d);
				dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, d);
				break;
			}
		}
	}
	if(!(wayMappingPos.ln2.endp.equals(ln2.getCoordinateN(wayMappingPos.ln2.endp_pos))))
	{
		double d = wayMappingPos.ln2.endp.distance(wayMappingPos.ln2.endp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}

	for(int i=wayMappingPos.ln2.endp_pos+1; i<ln2.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln2 = ln2.getCoordinateN(i);
		if(wayMappingPos.ln2.endp_mappingp_pos == ln1.getNumPoints()-1)
		{
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, pt_ln2.distance(ln1.getCoordinateN(ln1.getNumPoints()-1)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=wayMappingPos.ln2.endp_mappingp_pos+1; j<ln1.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
				d = std::min(d, sg_ln1.distance(pt_ln2));
			}
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, d);
		}
	}

	wayFeature.dis_average = dis_avg / num;
	wayFeature.dis_min = dis_min;
	wayFeature.dis_max = dis_max;
	wayFeature.dis_hausdorff_ln1 = dis_hausdorff_ln1;
	wayFeature.dis_hausdorff_ln2 = dis_hausdorff_ln2;
	return true;
}
