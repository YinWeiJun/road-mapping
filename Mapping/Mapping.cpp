/*
 * Mapping.cpp
 *
 *  Created on: 2018年5月4日
 *      Author: didi
 */

#include "way_feature_calculator.h"
#include <geos.h>
#include <cmath>
//#include <proj_api.h>

using namespace std;
const double ANGLE_THRESHOLD = M_PI/6;
const double LLBAND[6] = { 75, 60, 45, 30, 15, 0 };
const double LL2MC[][10] = { { -0.00157021024440, 1.113207020616939e+005,
		1.704480524535203e+015, -1.033898737604234e+016,
		2.611266785660388e+016, -3.514966917665370e+016,
		2.659570071840392e+016, -1.072501245418824e+016,
		1.800819912950474e+015, 82.50000000000000 }, { 8.277824516172526e-004,
		1.113207020463578e+005, 6.477955746671608e+008,
		-4.082003173641316e+009, 1.077490566351142e+010,
		-1.517187553151559e+010, 1.205306533862167e+010,
		-5.124939663577472e+009, 9.133119359512032e+008, 67.50000000000000 }, {
		0.00337398766765, 1.113207020202162e+005, 4.481351045890365e+006,
		-2.339375119931662e+007, 7.968221547186455e+007,
		-1.159649932797253e+008, 9.723671115602145e+007,
		-4.366194633752821e+007, 8.477230501135234e+006, 52.50000000000000 }, {
		0.00220636496208, 1.113207020209128e+005, 5.175186112841131e+004,
		3.796837749470245e+006, 9.920137397791013e+005,
		-1.221952217112870e+006, 1.340652697009075e+006,
		-6.209436990984312e+005, 1.444169293806241e+005, 37.50000000000000 }, {
		-3.441963504368392e-004, 1.113207020576856e+005,
		2.782353980772752e+002, 2.485758690035394e+006, 6.070750963243378e+003,
		5.482118345352118e+004, 9.540606633304236e+003,
		-2.710553267466450e+003, 1.405483844121726e+003, 22.50000000000000 },
		{-3.218135878613132e-004, 1.113207020701615e+005, 0.00369383431289,
				8.237256402795718e+005, 0.46104986909093,
				2.351343141331292e+003, 1.58060784298199, 8.77738589078284,
				0.37238884252424, 7.45000000000000 } };

WayFeatureCalculator::WayFeatureCalculator()
{
}

WayFeatureCalculator::~WayFeatureCalculator()
{
}

geos::geom::Coordinate WayFeatureCalculator::ll2mc(const geos::geom::Coordinate& coor)
{
	geos::geom::Coordinate temp;
	temp.x = coor.x;
	if(temp.x > 180.0) {
		temp.x = 180.0;
	} else if(temp.x < -180.0) {
		temp.x  = -180.0;
	}

	temp.y = coor.y;
	if (temp.y < 1E-7 && temp.y >= 0.0) {
		temp.y = 1E-7;
	} else if(temp.y < 0 && temp.y > -1.0E-7) {
		temp.y = -1E-7;
	} else if(temp.y > 74) {
		temp.y = 74;
	} else if(temp.y < -74) {
		temp.y = -74;
	}

	double factor[10] = { 0 };
	unsigned int i = 0;
	for (i = 0; i < sizeof(LLBAND) / sizeof(double); i++) {
		if (std::fabs(temp.y) > LLBAND[i]) {
			memcpy(factor, LL2MC[i], sizeof(factor));
			break;
		}
	}
	return _conv_(temp, factor);
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
			int in = GeosTool::pt2Line(tmp, sg_ln2.p0, sg_ln2.p1, prjc);
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
					while(j > 0)
					{
						geos::geom::LineSegment sg_ln2_tmp = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
						WayFeatureCalculator::calMappingPoint(startp, sg_ln1_tmp, sg_ln2_tmp, true);
						if(GeosTool::pt2Line(startp, sg_ln1_tmp.p0, sg_ln1_tmp.p1, unused) == 0)
							j -= 1;
						else
							break;
					}
					j += 1;
					sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
					WayFeatureCalculator::calMappingPoint(startp, sg_ln1_tmp, sg_ln2, true);
					if( GeosTool::pt2Line(startp, sg_ln1_tmp.p0, sg_ln1_tmp.p1, unused) == 0 && GeosTool::calculate_distance(sg_ln2.p0, startp) < dis_Threshold)
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
			if(GeosTool::pt2Line(startp, sg_ln1.p0, sg_ln1.p1, unused) == 0 && GeosTool::calculate_distance(sg_ln2.p0, startp) < dis_Threshold)
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
			int in = GeosTool::pt2Line(tmp, sg_ln2.p0, sg_ln2.p1, prjc);
			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
			double d_angle2 = M_PI;
			if (i < ln2.getNumPoints() - 1)
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
					while(j < ln2.getNumPoints())
					{
						geos::geom::LineSegment sg_ln2_tmp = geos::geom::LineSegment(ln2.getCoordinateN(j), ln2.getCoordinateN(j-1));
						WayFeatureCalculator::calMappingPoint(startp, sg_ln1_tmp, sg_ln2_tmp, true);
						if(GeosTool::pt2Line(startp, sg_ln1_tmp.p0, sg_ln1_tmp.p1, unused) == 0)
							j += 1;
						else
							break;
					}
					j -= 1;
					sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j), ln2.getCoordinateN(j-1));
					WayFeatureCalculator::calMappingPoint(startp, sg_ln1_tmp, sg_ln2, true);
					if(GeosTool::pt2Line(startp, sg_ln1_tmp.p0, sg_ln1_tmp.p1, unused) == 0 && GeosTool::calculate_distance(sg_ln2.p0, startp) < dis_Threshold)
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
			if(GeosTool::pt2Line(endp, sg_ln1.p0, sg_ln1.p1, unused) == 0 && GeosTool::calculate_distance(sg_ln2.p0, endp) < dis_Threshold)
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
bool WayFeatureCalculator::calHausdorffDis(WayHaustorffDis& wayHausdorff, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2)
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
			GeosTool::pt2Line(pt_ln1, sg_ln2.p0, sg_ln2.p1, unused);
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
			GeosTool::pt2Line(pt_ln2, sg_ln1.p0, sg_ln1.p1, unused);
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
			int in = GeosTool::pt2Line(pt_ln1, sg_ln2.p0, sg_ln2.p1, prjp);
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
			int in = GeosTool::pt2Line(pt_ln2, sg_ln1.p0, sg_ln1.p1, prjp);
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

	return true;
}

// ln1和ln2之间的豪斯多夫距离和匹配距离
bool WayFeatureCalculator::calMappingFeature(WayFeature& wayFeature, const geos::geom::LineString& ln1, const geos::geom::LineString& ln2, double dis_Threshold)
{
	geos::geom::Coordinate unused;
	WayMappingPos wayMappingPos;
	WayFeatureCalculator::calCoverOfLine(wayMappingPos, ln1, ln2, dis_Threshold);
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
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, GeosTool::calculate_distance(pt_ln1, ln2.getCoordinateN(0)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=1; j<=wayMappingPos.ln1.startp_mappingp_pos && j<ln1.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
				GeosTool::pt2Line(pt_ln1, sg_ln2.p0, sg_ln2.p1, unused);
				d = std::min(d,  GeosTool::calculate_distance(pt_ln1, unused));
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
		double d = GeosTool::calculate_distance(wayMappingPos.ln1.startp, wayMappingPos.ln1.startp_mappingp);
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
			int in = GeosTool::pt2Line(pt_ln1, sg_ln2.p0, sg_ln2.p1, prjp);
			if(in == 0)
			{
				num += 1;
				double d = GeosTool::calculate_distance(pt_ln1, prjp);
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
		double d = GeosTool::calculate_distance(wayMappingPos.ln1.endp, wayMappingPos.ln1.endp_mappingp);
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
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, GeosTool::calculate_distance(pt_ln1, ln2.getCoordinateN(ln2.getNumPoints()-1)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=wayMappingPos.ln1.endp_mappingp_pos+1; j<ln2.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
				GeosTool::pt2Line(pt_ln1, sg_ln2.p0, sg_ln2.p1, unused);
				d = std::min(d, GeosTool::calculate_distance(pt_ln1, unused));
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
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, GeosTool::calculate_distance(pt_ln2, ln1.getCoordinateN(0)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=1; j<=wayMappingPos.ln2.startp_mappingp_pos && j<ln2.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
				GeosTool::pt2Line(pt_ln2, sg_ln1.p0, sg_ln1.p1, unused);
				d = std::min(d, GeosTool::calculate_distance(pt_ln2, unused));
			}
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, d);
		}
	}

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
			int in = GeosTool::pt2Line(pt_ln2, sg_ln1.p0, sg_ln1.p1, prjp);
			if(in == 0)
			{
				num += 1;
				double d = GeosTool::calculate_distance(pt_ln2, prjp);
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
		double d = GeosTool::calculate_distance(wayMappingPos.ln2.endp, wayMappingPos.ln2.endp_mappingp);
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
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, GeosTool::calculate_distance(pt_ln2, ln1.getCoordinateN(ln1.getNumPoints()-1)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=wayMappingPos.ln2.endp_mappingp_pos+1; j<ln1.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
				GeosTool::pt2Line(pt_ln2, sg_ln1.p0, sg_ln1.p1, unused);
				d = std::min(d, GeosTool::calculate_distance(pt_ln2, unused));
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

// 坐标转换
geos::geom::Coordinate WayFeatureCalculator::_conv_(
		const geos::geom::Coordinate& fromPoint, double factor[])
{
	geos::geom::Coordinate toPoint;
	toPoint.x = factor[0] + factor[1] * fabs(fromPoint.x);
	double temp = fabs(fromPoint.y) / factor[9];
	toPoint.y = factor[2] + factor[3] * temp + factor[4] * temp * temp
			+ factor[5] * temp * temp * temp + factor[6] * temp * temp * temp
			* temp + factor[7] * temp * temp * temp * temp * temp + factor[8]
			* temp * temp * temp * temp * temp * temp;
	toPoint.x *= (fromPoint.x < 0 ? -1 : 1);
	toPoint.y *= (fromPoint.y < 0 ? -1 : 1);
	return toPoint;
}

