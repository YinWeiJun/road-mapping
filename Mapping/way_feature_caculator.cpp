/*
 * Mapping.cpp
 *
 *  Created on: 2018年5月4日
 *      Author: yinweijun
 */

#include "way_feature_calculator.h"
#include <geos.h>
#include <cmath>

using namespace std;
const double ANGLE_THRESHOLD = M_PI/3;

WayFeatureCalculator::WayFeatureCalculator()
{
}

WayFeatureCalculator::~WayFeatureCalculator()
{
}

// 求sg2 p0点作sg2的垂线，与sg1的交点 sg1与sg2的夹角小于30度
void WayFeatureCalculator::calMappingPoint(geos::geom::Coordinate& insp, const geos::geom::LineSegment& sg1, const geos::geom::LineSegment& sg2, bool first)
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
		b2 = pt.y - k2 * pt.x;
		if(sg1.p0.x == sg1.p1.x)
		{
			insp.x = sg1.p0.x;
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
bool WayFeatureCalculator::calMappingPointPos(WayMappingPt& wayMappingPt, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	geos::geom::Coordinate unused;
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
			int in = GeosTool::point_segment(tmp, sg_ln2.p0, sg_ln2.p1, prjc);
			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
			double d_angle2 = M_PI;
			if (i>1)
			{
				geos::geom::LineSegment sg_ln1_bef = geos::geom::LineSegment(ln1.getCoordinateN(i-2), ln1.getCoordinateN(i-1));
				double d_angle_tmp = fabs(sg_ln1_bef.angle() - angle_sg_ln2);
				if(d_angle_tmp > M_PI)
					d_angle_tmp = 2*M_PI - d_angle_tmp;
				d_angle2 = d_angle_tmp;
			}
			if(d_angle1 > M_PI)
				d_angle1 = 2*M_PI - d_angle1;
			if ( in == 0 && GeosTool::calculate_distance(tmp, prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < angle_Threshold))
			{
				wayMappingPt.startp = tmp;
				wayMappingPt.startp_pos = i-1;
				wayMappingPt.startp_mappingp = prjc;
				if (prjc.equals(sg_ln2.p0))
					wayMappingPt.startp_mappingp_pos = j - 1;
				else
					wayMappingPt.startp_mappingp_pos = j;
				if(1 < i && !prjc.equals(sg_ln2.p0))
				{
					geos::geom::LineSegment sg_ln1_tmp = geos::geom::LineSegment(ln1.getCoordinateN(i-2), ln1.getCoordinateN(i-1));
					geos::geom::Coordinate startp;
					j -= 1;
					while(j > 0)
					{
						geos::geom::LineSegment sg_ln2_tmp = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
						WayFeatureCalculator::calMappingPoint(startp, sg_ln1_tmp, sg_ln2_tmp, true);
						double d_angle_tmp = GeosTool::calculate_include_angle(sg_ln1_tmp.angle(), sg_ln2_tmp.angle());
						if(d_angle_tmp > M_PI)
							d_angle_tmp = 2*M_PI - d_angle_tmp;
						if(GeosTool::point_segment(startp, sg_ln1_tmp.p0, sg_ln1_tmp.p1, unused) == 0 && GeosTool::calculate_distance(startp, unused) < dis_Threshold && d_angle_tmp < angle_Threshold)
							j -= 1;
						else
							break;
					}
					j += 1;
					sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
					WayFeatureCalculator::calMappingPoint(startp, sg_ln1_tmp, sg_ln2, true);
					double d_angle = GeosTool::calculate_include_angle(sg_ln1_tmp.angle(), sg_ln2.angle());
					if(d_angle > M_PI)
						d_angle = 2*M_PI - d_angle;
					if( GeosTool::point_segment(startp, sg_ln1_tmp.p0, sg_ln1_tmp.p1, unused) == 0 && GeosTool::calculate_distance(sg_ln2.p0, startp) < dis_Threshold && d_angle < angle_Threshold)
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
		if(i == ln1.getNumPoints()-1 && !finded)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(0), ln2.getCoordinateN(1));
			geos::geom::Coordinate startp;
			WayFeatureCalculator::calMappingPoint(startp, sg_ln1, sg_ln2, true);
			double d_angle = fabs(sg_ln2.angle() - sg_ln1.angle());
			if(d_angle > M_PI)
				d_angle = 2*M_PI - d_angle;
			if(GeosTool::point_segment(startp, sg_ln1.p0, sg_ln1.p1, unused) == 0 && GeosTool::calculate_distance(sg_ln2.p0, startp) < dis_Threshold && d_angle < angle_Threshold)
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
			int in = GeosTool::point_segment(tmp, sg_ln2.p0, sg_ln2.p1, prjc);
			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
			double d_angle2 = M_PI;
			if (i < ln1.getNumPoints() - 1)
			{
				geos::geom::LineSegment sg_ln1_bef = geos::geom::LineSegment(ln1.getCoordinateN(i+1), ln1.getCoordinateN(i));
				double d_angle_tmp = fabs(sg_ln1_bef.angle() - angle_sg_ln2);
				if(d_angle_tmp > M_PI)
					d_angle_tmp = 2*M_PI - d_angle_tmp;
				d_angle2 = d_angle_tmp;
			}
			if(d_angle1 > M_PI)
				d_angle1 = 2*M_PI - d_angle1;
			if (in == 0 && GeosTool::calculate_distance(tmp, prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < angle_Threshold))
			{
				wayMappingPt.endp = tmp;
				wayMappingPt.endp_pos = i;
				wayMappingPt.endp_mappingp = prjc;
				if (prjc.equals(sg_ln2.p0))
					wayMappingPt.endp_mappingp_pos = j;
				else
					wayMappingPt.endp_mappingp_pos = j - 1;
				if(i < ln1.getNumPoints() - 1 && !prjc.equals(sg_ln2.p0))
				{
					geos::geom::LineSegment sg_ln1_tmp = geos::geom::LineSegment(ln1.getCoordinateN(i+1), ln1.getCoordinateN(i));
					geos::geom::Coordinate startp;
					j += 1;
					while(j < ln2.getNumPoints())
					{
						geos::geom::LineSegment sg_ln2_tmp = geos::geom::LineSegment(ln2.getCoordinateN(j), ln2.getCoordinateN(j-1));
						WayFeatureCalculator::calMappingPoint(startp, sg_ln1_tmp, sg_ln2_tmp, true);
						double d_angle_tmp = GeosTool::calculate_include_angle(sg_ln1_tmp.angle(), sg_ln2_tmp.angle());
						if (d_angle_tmp > M_PI)
							d_angle_tmp = 2*M_PI - d_angle_tmp;
						if(GeosTool::point_segment(startp, sg_ln1_tmp.p0, sg_ln1_tmp.p1, unused) == 0 && GeosTool::calculate_distance(startp, unused) < dis_Threshold && d_angle_tmp < angle_Threshold)
							j += 1;
						else
							break;
					}
					j -= 1;
					sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j), ln2.getCoordinateN(j-1));
					WayFeatureCalculator::calMappingPoint(startp, sg_ln1_tmp, sg_ln2, true);
					double d_angle = GeosTool::calculate_include_angle(sg_ln1_tmp.angle(), sg_ln2.angle());
					if (d_angle > M_PI)
						d_angle = 2*M_PI - d_angle;
					if(GeosTool::point_segment(startp, sg_ln1_tmp.p0, sg_ln1_tmp.p1, unused) == 0 && GeosTool::calculate_distance(sg_ln2.p0, startp) < dis_Threshold && d_angle < angle_Threshold)
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
		if(1 == i && finded2)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(ln2.getNumPoints()-1), ln2.getCoordinateN(ln2.getNumPoints()-2));
			geos::geom::Coordinate endp;
			WayFeatureCalculator::calMappingPoint(endp, sg_ln1, sg_ln2, dis_Threshold);
			double d_angle = fabs(sg_ln2.angle() - sg_ln1.angle());
			if(d_angle > M_PI)
				d_angle = 2*M_PI - d_angle;
			if(GeosTool::point_segment(endp, sg_ln1.p0, sg_ln1.p1, unused) == 0 && GeosTool::calculate_distance(sg_ln2.p0, endp) < dis_Threshold && d_angle < angle_Threshold)
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

	return finded && (!finded2);
}

// ln1与ln2, 相互到第一个匹配点和最后一个匹配点，以及对应到映射点
bool WayFeatureCalculator::calCoverOfLine(WayMappingPos& wayMappingPos, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	WayMappingPt wayMappingPt_ln1;
	bool mapping_ln1 = WayFeatureCalculator::calMappingPointPos(wayMappingPt_ln1, ln1, ln2, dis_Threshold);
	wayMappingPos.ln1 = wayMappingPt_ln1;
	wayMappingPos.mapping_ln1 = mapping_ln1;
	WayMappingPt wayMappingPt_ln2;
	bool mapping_ln2 = WayFeatureCalculator::calMappingPointPos(wayMappingPt_ln2, ln2, ln1, dis_Threshold);
	wayMappingPos.ln2 = wayMappingPt_ln2;
	wayMappingPos.mapping_ln2 = mapping_ln2;

	if(!mapping_ln1 && mapping_ln2)
	{
		wayMappingPos.ln1.startp = wayMappingPos.ln2.startp_mappingp;
		wayMappingPos.ln1.startp_pos = wayMappingPos.ln2.startp_mappingp_pos;
		wayMappingPos.ln1.startp_mappingp = wayMappingPos.ln2.startp;
		wayMappingPos.ln1.startp_mappingp_pos = wayMappingPos.ln2.startp_pos;
		wayMappingPos.ln1.endp = wayMappingPos.ln2.endp_mappingp;
		wayMappingPos.ln1.endp_pos = wayMappingPos.ln2.endp_mappingp_pos;
		wayMappingPos.ln1.endp_mappingp = wayMappingPos.ln2.endp;
		wayMappingPos.ln1.endp_mappingp_pos = wayMappingPos.ln2.endp_pos;
	}
	else if(mapping_ln1 && !mapping_ln2)
	{
		wayMappingPos.ln2.startp = wayMappingPos.ln1.startp_mappingp;
		wayMappingPos.ln2.startp_pos = wayMappingPos.ln1.startp_mappingp_pos;
		wayMappingPos.ln2.startp_mappingp = wayMappingPos.ln1.startp;
		wayMappingPos.ln2.startp_mappingp_pos = wayMappingPos.ln1.startp_pos;
		wayMappingPos.ln2.endp = wayMappingPos.ln1.endp_mappingp;
		wayMappingPos.ln2.endp_pos = wayMappingPos.ln1.endp_mappingp_pos;
		wayMappingPos.ln2.endp_mappingp = wayMappingPos.ln1.endp;
		wayMappingPos.ln2.endp_mappingp_pos = wayMappingPos.ln1.endp_pos;
	}

//	std::cerr << "\n==========================================" << std::endl;
//	std::cerr << "	Line1:	" + ln1.toString() << std::endl;
//	std::cerr << "	Line2:	" + ln2.toString() << std::endl;
//	std::cerr << "	Length1:	" + std::to_string(ln1.getLength()) << std::endl;
//	std::cerr << "	Length2:	" + std::to_string(ln2.getLength()) << std::endl;
//	std::cerr << "----------------------------------" << std::endl;
//	std::cerr << "	StartP:	" + wayMappingPos.ln1.startp.toString() + "	Pos:	" + std::to_string(wayMappingPos.ln1.startp_pos) + "	StartMappingP:	" + wayMappingPos.ln1.startp_mappingp.toString() + "	Pos:	" + std::to_string(wayMappingPos.ln1.startp_mappingp_pos) << std::endl;
//	std::cerr << "----------------------------------" << std::endl;
//	std::cerr << "	EndP:	" + wayMappingPos.ln1.endp.toString() + "	Pos:	" + std::to_string(wayMappingPos.ln1.endp_pos) + "	EndMappingP:	" + wayMappingPos.ln1.endp_mappingp.toString() + "	Pos:	" + std::to_string(wayMappingPos.ln1.endp_mappingp_pos) << std::endl;
//	std::cerr << "----------------------------------" << std::endl;
//	std::cerr << "	StartP:	" + wayMappingPos.ln2.startp.toString() + "	Pos:	" + std::to_string(wayMappingPos.ln2.startp_pos) + "	StartMappingP:	" + wayMappingPos.ln2.startp_mappingp.toString() + "	Pos:	" + std::to_string(wayMappingPos.ln2.startp_mappingp_pos) << std::endl;
//	std::cerr << "----------------------------------" << std::endl;
//	std::cerr << "	EndP:	" + wayMappingPos.ln2.endp.toString() + "	Pos:	" + std::to_string(wayMappingPos.ln2.endp_pos) + "	EndMappingP:	" + wayMappingPos.ln2.endp_mappingp.toString() + "	Pos:	" + std::to_string(wayMappingPos.ln2.endp_mappingp_pos) << std::endl;
//	std::cerr << "----------------------------------" << std::endl;

	if(!mapping_ln1 && !mapping_ln2)
		return false;
	else
		return true;
}

// ln1与ln2之间点豪斯多夫距离 包括ln1到ln2和ln2到ln1的两个豪斯多夫距离的最大和最小值
bool WayFeatureCalculator::calHausdorffDis(WayHausdorffDis& wayHausdorff, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2)
{
	geos::geom::Coordinate unused;
	double dis_hausdorff_ln1 = 0;
	int j = 1;
	for(int i=0; i<ln1.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		double d = GeosTool::calculate_distance(ln2.getCoordinateN(j-1),pt_ln1);
		for( ; j<ln2.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
			GeosTool::point_segment(pt_ln1, sg_ln2.p0, sg_ln2.p1, unused);
			double dis = GeosTool::calculate_distance(pt_ln1, unused);
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
		double d = GeosTool::calculate_distance(ln1.getCoordinateN(j-1), pt_ln2);
		for( ; j<ln1.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
			GeosTool::point_segment(pt_ln2, sg_ln1.p0, sg_ln1.p1, unused);
			double dis = GeosTool::calculate_distance(pt_ln2, unused);
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

// 计算ln1与ln2的豪斯多夫角度
bool WayFeatureCalculator::calHausdorffDir(WayHausdorffDir& wayHausdorffDir, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2)
{
	double dir_hausdorff_ln1(0.);

	for(int i=1; i<ln1.getNumPoints(); i++)
	{
		geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(i-1), ln1.getCoordinateN(i));
		double angle_ln1 = sg_ln1.angle();
		double dir(DBL_MAX);
		for(int j = 1; j<ln2.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
			double angle_ln2 = sg_ln2.angle();
			dir = std::min(dir, GeosTool::calculate_include_angle(angle_ln1, angle_ln2));
		}
		dir_hausdorff_ln1 = std::max(dir_hausdorff_ln1, dir);
	}

	double dir_hausdorff_ln2(0.);

	for(int i=1; i<ln2.getNumPoints(); i++)
	{
		geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(i-1), ln2.getCoordinateN(i));
		double angle_ln2 = sg_ln2.angle();
		double dir(DBL_MAX);
		for(int j = 1; j<ln1.getNumPoints(); j++)
		{
			geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
			double angle_ln1 = sg_ln1.angle();
			dir = std::min(dir, GeosTool::calculate_include_angle(angle_ln2, angle_ln1));
		}
		dir_hausdorff_ln2 = std::max(dir_hausdorff_ln2, dir);
	}
	wayHausdorffDir.dir_hausdorff_ln1 = dir_hausdorff_ln1;
	wayHausdorffDir.dir_hausdorff_ln2 = dir_hausdorff_ln2;
	return true;
}

// ln1与ln2匹配点间的距离 包括最大／最小／平均
bool WayFeatureCalculator::calMappingDis(WayMappingDis& wayMappingDis, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	geos::geom::Coordinate unused;
	WayMappingPos wayMappingPos;
	WayFeatureCalculator::calCoverOfLine(wayMappingPos, ln1, ln2, dis_Threshold);
	int num = 0;
	double dis_avg = 0;
	double dis_min = DBL_MAX;
	double dis_max = 0;
	int j = wayMappingPos.ln1.startp_mappingp_pos;
	if(0 == j)
		j += 1;
	if(!(wayMappingPos.ln1.startp.equals(ln1.getCoordinateN(wayMappingPos.ln1.startp_pos))))
	{
		double d = GeosTool::calculate_distance(wayMappingPos.ln1.startp, wayMappingPos.ln1.startp_mappingp);
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
			int in = GeosTool::point_segment(pt_ln1, sg_ln2.p0, sg_ln2.p1, prjp);
			if(in == 0)
			{
				num += 1;
				double d = GeosTool::calculate_distance(pt_ln1, prjp);
				dis_avg += d;
				dis_min = std::min(dis_min, d);
				dis_max = std::max(dis_max, d);
				break;
			}
		}
	}
	if(!(wayMappingPos.ln1.endp.equals(ln1.getCoordinateN(wayMappingPos.ln1.endp_pos))))
	{
		double d = GeosTool::calculate_distance(wayMappingPos.ln1.endp, wayMappingPos.ln1.endp_mappingp);
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
		double d = GeosTool::calculate_distance(wayMappingPos.ln2.startp, wayMappingPos.ln2.startp_mappingp);
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
			int in = GeosTool::point_segment(pt_ln2, sg_ln1.p0, sg_ln1.p1, prjp);
			if(in == 0)
			{
				num += 1;
				double d = GeosTool::calculate_distance(pt_ln2, prjp);
				dis_avg += d;
				dis_min = std::min(dis_min, d);
				dis_max = std::max(dis_max, d);
				break;
			}
		}
	}
	if(!(wayMappingPos.ln2.endp.equals(ln2.getCoordinateN(wayMappingPos.ln2.endp_pos))))
	{
		double d = GeosTool::calculate_distance(wayMappingPos.ln2.endp, wayMappingPos.ln2.endp_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}

	wayMappingDis.dis_avg = dis_avg / num;
	wayMappingDis.dis_min = dis_min;
	wayMappingDis.dis_max = dis_max;
	wayMappingDis.dis_mapping_sta = GeosTool::calculate_distance(wayMappingPos.ln1.startp, wayMappingPos.ln1.startp_mappingp);
	wayMappingDis.dis_mapping_end = GeosTool::calculate_distance(wayMappingPos.ln1.endp, wayMappingPos.ln1.endp_mappingp);
	return true;
}

// 计算ln1和ln2的匹配长度
bool WayFeatureCalculator::calMappingLength(WayCoverLength& wayCoverLength,
		const geos::geom::LineString& ln1, const geos::geom::LineString& ln2,
		double dis_Threshold)
{
	WayMappingPos wayMappingPos;
	if (!WayFeatureCalculator::calCoverOfLine(wayMappingPos, ln1, ln2, dis_Threshold))
		return false;

	double len_ln1_cover = 0.;
	int left = wayMappingPos.ln1.startp_pos;
	int right = wayMappingPos.ln1.endp_pos;
	if(left != -1 && right != -1)
	{
		if(!wayMappingPos.ln1.startp.equals(ln1.getCoordinateN(wayMappingPos.ln1.startp_pos)) && wayMappingPos.ln1.startp_pos <= wayMappingPos.ln1.endp_pos)
		{
			len_ln1_cover += GeosTool::calculate_distance(wayMappingPos.ln1.startp, ln1.getCoordinateN(wayMappingPos.ln1.startp_pos));
		}
		while(left < right)
		{

			len_ln1_cover += GeosTool::calculate_distance(ln1.getCoordinateN(left), ln1.getCoordinateN(left+1));
			left += 1;
		}
		if(!wayMappingPos.ln1.endp.equals(ln1.getCoordinateN(wayMappingPos.ln1.endp_pos)) && wayMappingPos.ln1.startp_pos <= wayMappingPos.ln1.endp_pos)
		{
			len_ln1_cover += GeosTool::calculate_distance(wayMappingPos.ln1.endp, ln1.getCoordinateN(wayMappingPos.ln1.endp_pos));
		}
		if(left > right)
		{
			len_ln1_cover += GeosTool::calculate_distance(wayMappingPos.ln1.startp, wayMappingPos.ln1.endp);
		}
	}

	double len_ln2_cover = 0.;
	left = wayMappingPos.ln2.startp_pos;
	right = wayMappingPos.ln2.endp_pos;
	if(left != -1 && right != -1)
	{
		if(!wayMappingPos.ln2.startp.equals(ln2.getCoordinateN(wayMappingPos.ln2.startp_pos)) && wayMappingPos.ln2.startp_pos <= wayMappingPos.ln2.endp_pos)
		{
			len_ln2_cover += GeosTool::calculate_distance(wayMappingPos.ln2.startp, ln2.getCoordinateN(wayMappingPos.ln2.startp_pos));
		}
		while(left < right)
		{
			len_ln2_cover += GeosTool::calculate_distance(ln2.getCoordinateN(left), ln2.getCoordinateN(left+1));
			left += 1;
		}
		if(!wayMappingPos.ln2.endp.equals(ln2.getCoordinateN(wayMappingPos.ln2.endp_pos)) && wayMappingPos.ln2.startp_pos <= wayMappingPos.ln2.endp_pos)
		{
			len_ln2_cover += GeosTool::calculate_distance(wayMappingPos.ln2.endp, ln2.getCoordinateN(wayMappingPos.ln2.endp_pos));
		}
		if(left > right)
		{
			len_ln2_cover += GeosTool::calculate_distance(wayMappingPos.ln2.startp, wayMappingPos.ln2.endp);
		}
	}

	wayCoverLength.len_ln1_cover = len_ln1_cover;
	wayCoverLength.len_ln2_cover = len_ln2_cover;

	if(wayCoverLength.len_ln1_cover < 1e-7 || wayCoverLength.len_ln2_cover < 1e-7)
	{
		wayCoverLength.len_ln1_cover += wayCoverLength.len_ln2_cover;
		wayCoverLength.len_ln2_cover = wayCoverLength.len_ln1_cover;
	}

	double angle_diff_cover(DBL_MAX);
	geos::geom::LineSegment sg1 = geos::geom::LineSegment(wayMappingPos.ln1.startp, wayMappingPos.ln1.endp);
	geos::geom::LineSegment sg2 = geos::geom::LineSegment(wayMappingPos.ln2.startp, wayMappingPos.ln2.endp);
	angle_diff_cover = GeosTool::calculate_include_angle(sg1.angle(), sg2.angle());
	wayCoverLength.angle_diff_cover = angle_diff_cover;


//	std::cerr << angle_diff_cover << "," << wayMappingPos.ln1.startp << "," << wayMappingPos.ln1.endp << "," << wayMappingPos.ln2.startp << "," << wayMappingPos.ln2.endp << std::endl;
	return true;
}

// ln1和ln2之间的多个特征
bool WayFeatureCalculator::calMappingFeature(WayFeature& wayFeature, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	WayCoverLength wayCoverLength;
	bool ret_wfc = WayFeatureCalculator::calMappingLength(wayCoverLength, ln1, ln2, dis_Threshold);
	if(ret_wfc)
	{
		wayFeature.length_mapping_ln1 = wayCoverLength.len_ln1_cover;
		wayFeature.length_mapping_ln2 = wayCoverLength.len_ln2_cover;
		wayFeature.angle_diff_mapping = wayCoverLength.angle_diff_cover;
	}
	WayHausdorffDis wayHausdorffDis;
	bool ret_wfc2 = WayFeatureCalculator::calHausdorffDis(wayHausdorffDis, ln1, ln2);
	if(ret_wfc2)
	{
		wayFeature.dis_hausdorff_ln1 = wayHausdorffDis.dis_hausdorff_ln1;
		wayFeature.dis_hausdorff_ln2 = wayHausdorffDis.dis_hausdorff_ln2;
	}
	WayHausdorffDir wayHausdorffDir;
	bool ret_wfc3 = WayFeatureCalculator::calHausdorffDir(wayHausdorffDir, ln1, ln2);
	if(ret_wfc3)
	{
		wayFeature.angle_hausdorff_ln1 = wayHausdorffDir.dir_hausdorff_ln1;
		wayFeature.angle_hausdorff_ln2 = wayHausdorffDir.dir_hausdorff_ln2;
	}
	WayMappingDis wayMappingDis;
	bool ret_wfc4 = WayFeatureCalculator::calMappingDis(wayMappingDis, ln1, ln2, dis_Threshold);
	if(ret_wfc4)
	{
		wayFeature.dis_average = wayMappingDis.dis_avg;
		wayFeature.dis_min = wayMappingDis.dis_min;
		wayFeature.dis_max = wayMappingDis.dis_max;
		wayFeature.dis_mapping_sta = wayMappingDis.dis_mapping_sta;
		wayFeature.dis_mapping_end = wayMappingDis.dis_mapping_end;
	}
	return ret_wfc && ret_wfc2 && ret_wfc3 && ret_wfc4;
}

