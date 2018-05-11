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

Mapping::Mapping()
{
}

Mapping::~Mapping()
{
}

geos::geom::Coordinate Mapping::calPt_pt2Segment(geos::geom::Point& p, geos::geom::Point& p1, geos::geom::Point& p2)
{
	geos::geom::Coordinate vector10(p.getX() - p1.getX(), p.getY() - p1.getY());
	geos::geom::Coordinate vector12(p2.getX() - p1.getX(), p2.getY() - p1.getY());

	double divisor = vector10.x * vector12.x + vector10.y * vector12.y;
	double Divisor = vector12.x * vector12.x + vector12.y * vector12.y;
	double r = divisor / Divisor;

	geos::geom::Coordinate resc;
	if (r <= 0)
	{
		resc.x = p1.getX();
		resc.y = p1.getY();
	}
	else if (r >= 1)
	{
		resc.x = p2.getX();
		resc.y = p2.getY();
	}
	else
	{
		resc.x = (1 - r) * p1.getX() + r * p2.getX();
		resc.y = (1 - r) * p1.getY() + r * p2.getY();
	}
	return resc;
}


geos::geom::Coordinate Mapping::calPt_pt2Segment(geos::geom::Point* p, geos::geom::Point* p1, geos::geom::Point* p2)
{
	geos::geom::Coordinate vector10(p->getX() - p1->getX(), p->getY() - p1->getY());
	geos::geom::Coordinate vector12(p2->getX() - p1->getX(), p2->getY() - p1->getY());

	double divisor = vector10.x * vector12.x + vector10.y * vector12.y;
	double Divisor = vector12.x * vector12.x + vector12.y * vector12.y;
	double r = divisor / Divisor;

	geos::geom::Coordinate resc;
	if (r <= 0)
	{
		resc.x = p1->getX();
		resc.y = p1->getY();
	}
	else if (r >= 1)
	{
		resc.x = p2->getX();
		resc.y = p2->getY();
	}
	else
	{
		resc.x = (1 - r) * p1->getX() + r * p2->getX();
		resc.y = (1 - r) * p1->getY() + r * p2->getY();
	}
	return resc;
}


geos::geom::Coordinate Mapping::calPt_Pt2Line(geos::geom::Point& p,
		geos::geom::LineString& ls) {
	geos::geom::Coordinate resc;
		double dis = 0;
		geos::geom::GeometryFactory::unique_ptr factory = geos::geom::GeometryFactory::create();
		for(int i=1; i<ls.getNumPoints(); i++)
		{
			geos::geom::CoordinateArraySequence *cas = new geos::geom::CoordinateArraySequence();
			cas->add(ls.getCoordinateN(i-1));
			cas->add(ls.getCoordinateN(i));
			geos::geom::LineString *l = factory->createLineString(cas);
			double d = p.distance(l);
			if (0 == dis)
			{
				dis = d;
				resc = Mapping::calPt_pt2Segment(p, *(ls.getPointN(i-1)), *(ls.getPointN(i)));
			}
			else if (d < dis)
			{
				dis = d;
				resc = Mapping::calPt_pt2Segment(p, *(ls.getPointN(i-1)), *(ls.getPointN(i)));
			}
		}
		return resc;
}


geos::geom::Coordinate Mapping::calPt_Pt2Line(geos::geom::Point* p,
		geos::geom::LineString* ls) {
	geos::geom::Coordinate resc;
	double dis = 0;
	geos::geom::GeometryFactory::unique_ptr factory = geos::geom::GeometryFactory::create();
	for(int i=1; i<ls->getNumPoints(); i++)
	{
		geos::geom::CoordinateArraySequence *cas = new geos::geom::CoordinateArraySequence();
		cas->add(ls->getCoordinateN(i-1));
		cas->add(ls->getCoordinateN(i));
		geos::geom::LineString *l = factory->createLineString(cas);
		double d = p->distance(l);
		if (0 == dis)
		{
			dis = d;
			resc = Mapping::calPt_pt2Segment(p, ls->getPointN(i-1), ls->getPointN(i));
		}
		else if (d < dis)
		{
			dis = d;
			resc = Mapping::calPt_pt2Segment(p, ls->getPointN(i-1), ls->getPointN(i));
		}
	}
	return resc;
}

// 求sg2 p0点作sg2的垂线，与sg1的交点 sg1与sg2的夹角小于30度
geos::geom::Coordinate Mapping::calMappingPoint(geos::geom::LineSegment& sg1, geos::geom::LineSegment& sg2)
{
	geos::geom::Coordinate insp = geos::geom::Coordinate();

	double k1 = 0;
	double b1 = 0;
	double k2 = 0;
	double b2 = 0;
	if(sg2.p0.y == sg2.p1.y)
	{
		insp.x = sg2.p0.x;
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
		insp.y = sg2.p0.y;
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

	return insp;
}

// ln1匹配到ln2的第一个点坐标及其位置，以及映射到了ln2上到坐标和位置
void Mapping::calMappingStartPoint(geos::geom::LineString& ln1, geos::geom::LineString& ln2, double dis_Threshold, CoverInfo& ci)
{
	double angle_Threshold = M_PI/6;

	int i = 1;
	for( ; i<ln1.getNumPoints() && ci.covered; i++)
	{
		geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(i-1), ln1.getCoordinateN(i));
		geos::geom::Coordinate tmp = ln1.getCoordinateN(i-1);
		double angle_sg_ln1 = sg_ln1.angle();
		int j = 1;
		for( ; j<ln2.getNumPoints(); j++)
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
			if (sg_ln2.distance(prjc) < 1e-7 && tmp.distance(prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < d_angle2 < angle_Threshold))
			{
				if(1 < i and !prjc.equals(sg_ln2.p0))
				{
					geos::geom::LineSegment sg_ln1_tmp = geos::geom::LineSegment(ln1.getCoordinateN(i-2), ln1.getCoordinateN(i-1));
					geos::geom::Coordinate startp = Mapping::calMappingPoint(sg_ln1_tmp, sg_ln2);
					if(sg_ln1_tmp.distance(startp) < 1e-7 && sg_ln2.p0.distance(startp) < dis_Threshold)
					{
						ci.startp_ln1 = startp;
						ci.startp_ln1_pos = i-1;
						ci.startp_ln1_mappingp = sg_ln2.p0;
						ci.startp_ln1_mappingp_pos = j-1;
					}
				}
				else
				{
					ci.startp_ln1 = tmp;
					ci.startp_ln1_pos = i-1;
					ci.startp_ln1_mappingp = prjc;
					if (prjc.equals(sg_ln2.p0))
						ci.startp_ln1_mappingp_pos = j - 1;
					else
						ci.startp_ln1_mappingp_pos = j;
				}
				j = ln2.getNumPoints()+1;
				i = ln1.getNumPoints()+1;
			}
		}
	}
	if(i==ln1.getNumPoints())
	{
		ci.covered = false;
	}
}

// ln1与ln2, 相互到第一个匹配点和最后一个匹配点，以及对应到映射点
CoverInfo Mapping::calCoverOfLine(geos::geom::LineString& ln1, geos::geom::LineString& ln2, double dis_Threshold)
{
	double angle_Threshold = M_PI/6;

	CoverInfo coverInfo = CoverInfo();
	CoverInfo ci_start_ln1 = CoverInfo();
	Mapping::calMappingStartPoint(ln1, ln2, dis_Threshold, ci_start_ln1);
	coverInfo.startp_ln1 = ci_start_ln1.startp_ln1;
	coverInfo.startp_ln1_pos = ci_start_ln1.startp_ln1_pos;
	coverInfo.startp_ln1_mappingp = ci_start_ln1.startp_ln1_mappingp;
	coverInfo.startp_ln1_mappingp_pos = ci_start_ln1.startp_ln1_mappingp_pos;

	CoverInfo ci_start_ln2 = CoverInfo();
	Mapping::calMappingStartPoint(ln2, ln1, dis_Threshold, ci_start_ln2);
	coverInfo.startp_ln2 = ci_start_ln2.startp_ln1;
	coverInfo.startp_ln2_pos = ci_start_ln2.startp_ln1_pos;
	coverInfo.startp_ln2_mappingp = ci_start_ln2.startp_ln1_mappingp;
	coverInfo.startp_ln2_mappingp_pos = ci_start_ln2.startp_ln1_mappingp_pos;

	geos::geom::Geometry* geo_ln1 = ln1.reverse();
	geos::geom::LineString* ln1_reverse = dynamic_cast<geos::geom::LineString*>(geo_ln1);

	geos::geom::Geometry* geo_ln2 = ln2.reverse();
	geos::geom::LineString* ln2_reverse = dynamic_cast<geos::geom::LineString*>(geo_ln2);

	CoverInfo ci_end_ln1 = CoverInfo();
	Mapping::calMappingStartPoint(*ln1_reverse, *ln2_reverse, dis_Threshold, ci_end_ln1);
	coverInfo.endp_ln1 = ci_end_ln1.startp_ln1;
	coverInfo.endp_ln1_pos = ln1.getNumPoints() - 1 - ci_end_ln1.startp_ln1_pos;
	coverInfo.endp_ln1_mappingp = ci_end_ln1.startp_ln1_mappingp;
	coverInfo.endp_ln1_mappingp_pos = ln2.getNumPoints() - 1 - ci_end_ln1.startp_ln1_mappingp_pos;

	CoverInfo ci_end_ln2 = CoverInfo();
	Mapping::calMappingStartPoint(*ln2_reverse, *ln1_reverse, dis_Threshold, ci_end_ln2);
	coverInfo.endp_ln2 = ci_end_ln2.startp_ln1;
	coverInfo.endp_ln2_pos = ln1.getNumPoints() - 1 - ci_end_ln2.startp_ln1_pos;
	coverInfo.endp_ln2_mappingp = ci_end_ln2.startp_ln1_mappingp;
	coverInfo.endp_ln2_mappingp_pos = ln2.getNumPoints() - 1 - ci_end_ln2.startp_ln1_mappingp_pos;

	return coverInfo;

//	CoverInfo ci;
//	int i = 1;
//	for( ; i<ln1.getNumPoints() && ci.covered; i++)
//	{
//		geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(i-1), ln1.getCoordinateN(i));
//		geos::geom::Coordinate tmp = ln1.getCoordinateN(i-1);
//		double angle_sg_ln1 = sg_ln1.angle();
//		int j = 1;
//		for( ; j<ln2.getNumPoints(); j++)
//		{
//			geos::geom::Coordinate prjc = geos::geom::Coordinate();
//			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
//			double angle_sg_ln2 = sg_ln2.angle();
//			sg_ln2.project(tmp, prjc);
//			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
//			double d_angle2 = 0;
//			if (i>1)
//			{
//				geos::geom::LineSegment sg_ln1_bef = geos::geom::LineSegment(ln1.getCoordinateN(i-2), ln1.getCoordinateN(i-1));
//				d_angle2 = fabs(sg_ln1_bef.angle() - angle_sg_ln2);
//			}
//			if (sg_ln2.distance(prjc) < 1e-7 && tmp.distance(prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < d_angle2 < angle_Threshold))
//			{
//				if(1==i || prjc.equals(sg_ln2.p0))
//				{
//					ci.startp_ln1 = tmp;
//					ci.startp_ln1_pos = i-1;
//					ci.startp_ln1_mappingp = prjc;
//					ci.startp_ln1_mappingp_pos = j;
//				}
//				else
//				{
//					geos::geom::LineSegment sg_ln1_tmp = geos::geom::LineSegment(ln1.getCoordinateN(i-2), ln1.getCoordinateN(i-1));
//					geos::geom::Coordinate startp = Mapping::calMappingPoint(sg_ln1_tmp, sg_ln2);
//					if(sg_ln1_tmp.distance(startp) < 1e-7 && sg_ln2.p0.distance(startp) < dis_Threshold)
//					{
//						ci.startp_ln1 = startp;
//						ci.startp_ln1_pos = i-1;
//						ci.startp_ln1_mappingp = sg_ln2.p0;
//						ci.startp_ln1_mappingp_pos = j-1;
//					}
//					else
//					{
//						ci.startp_ln1 = tmp;
//						ci.startp_ln1_pos = i-1;
//						ci.startp_ln1_mappingp = prjc;
//						ci.startp_ln1_mappingp_pos = j;
//					}
//				}
//				j = ln2.getNumPoints()+1;
//				i = ln1.getNumPoints()+1;
//			}
//		}
//	}
//	if(i==ln1.getNumPoints())
//	{
//		ci.covered = false;
//	}
//
//	i = ln1.getNumPoints()-1;
//	for( ; i>0 && ci.covered; i--)
//	{
//		geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(i), ln1.getCoordinateN(i-1));
//		geos::geom::Coordinate tmp = ln1.getCoordinateN(i);
//		double angle_sg_ln1 = sg_ln1.angle();
//		int j = ln2.getNumPoints() - 1;
//		for( ; j>0; j--)
//		{
//			geos::geom::Coordinate prjc = geos::geom::Coordinate();
//			geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j), ln2.getCoordinateN(j-1));
//			double angle_sg_ln2 = sg_ln2.angle();
//			sg_ln2.project(tmp, prjc);
//			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
//			double d_angle2 = 0;
//			if (i < ln1.getNumPoints() - 1)
//			{
//				geos::geom::LineSegment sg_ln1_bef = geos::geom::LineSegment(ln1.getCoordinateN(i+1), ln1.getCoordinateN(i));
//				d_angle2 = fabs(sg_ln1_bef.angle() - angle_sg_ln2);
//			}
//			if (sg_ln2.distance(prjc) < 1e-7  && tmp.distance(prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < angle_Threshold))
//			{
//				if(i == ln1.getNumPoints()-1 or prjc.equals(sg_ln2.p0))
//				{
//					ci.endp_ln1 = tmp;
//					ci.endp_ln1_pos = i;
//					ci.endp_ln1_mappingp = prjc;
//					ci.endp_ln1_mappingp_pos = j-1;
//				}
//				else
//				{
//					geos::geom::LineSegment sg_ln1_tmp = geos::geom::LineSegment(ln1.getCoordinateN(i+1), ln1.getCoordinateN(i));
//					geos::geom::Coordinate endp = Mapping::calMappingPoint(sg_ln1_tmp ,sg_ln2);
//					if(sg_ln1_tmp.distance(endp) < 1e-7 && sg_ln2.p0.distance(endp) < dis_Threshold)
//					{
//						ci.endp_ln1 = endp;
//						ci.endp_ln1_pos = i;
//						ci.endp_ln1_mappingp = sg_ln2.p0;
//						ci.endp_ln1_mappingp_pos = j;
//					}
//					else
//					{
//						ci.endp_ln1 = tmp;
//						ci.endp_ln1_pos = i;
//						ci.endp_ln1_mappingp = prjc;
//						ci.endp_ln1_mappingp_pos = j-1;
//					}
//				}
//				j = -1;
//				i = -1;
//			}
//		}
//	}
//	if(i==0)
//	{
//		ci.covered = false;
//	}
//
//	i = 1;
//	for( ; i<ln2.getNumPoints() && ci.covered; i++)
//	{
//		geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(i-1), ln2.getCoordinateN(i));
//		geos::geom::Coordinate tmp = ln2.getCoordinateN(i-1);
//		double angle_sg_ln2 = sg_ln2.angle();
//		int j = 1;
//		for( ; j<ln1.getNumPoints(); j++)
//		{
//			geos::geom::Coordinate prjc = geos::geom::Coordinate();
//			geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
//			double angle_sg_ln1 = sg_ln1.angle();
//			sg_ln1.project(tmp, prjc);
//			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
//			double d_angle2 = 0;
//			if (i>1)
//			{
//				geos::geom::LineSegment sg_ln2_bef = geos::geom::LineSegment(ln2.getCoordinateN(i-2), ln2.getCoordinateN(i-1));
//				d_angle2 = fabs(sg_ln2_bef.angle() - angle_sg_ln1);
//			}
//			if (sg_ln1.distance(prjc) < 1e-7  && tmp.distance(prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < angle_Threshold))
//			{
//				if(1 == i || prjc.equals(sg_ln1.p0))
//				{
//					ci.startp_ln2 = tmp;
//					ci.startp_ln2_pos = i-1;
//					ci.startp_ln2_mappingp = prjc;
//					ci.startp_ln2_mappingp_pos = j;
//				}
//				else
//				{
//					geos::geom::LineSegment sg_ln2_tmp = geos::geom::LineSegment(ln2.getCoordinateN(i-2), ln2.getCoordinateN(i-1));
//					geos::geom::Coordinate startp = Mapping::calMappingPoint(sg_ln2_tmp, sg_ln1);
//					if(sg_ln2_tmp.distance(startp) < 1e-7 && sg_ln2.p0.distance(startp) < dis_Threshold)
//					{
//						ci.startp_ln2 = startp;
//						ci.startp_ln2_pos = i-1;
//						ci.startp_ln2_mappingp = sg_ln1.p0;
//						ci.startp_ln2_mappingp_pos = j-1;
//					}
//					else
//					{
//						ci.startp_ln2 = tmp;
//						ci.startp_ln2_pos = i-1;
//						ci.startp_ln2_mappingp = prjc;
//						ci.startp_ln2_mappingp_pos = j;
//					}
//				}
//				j = ln1.getNumPoints()+1;
//				i = ln2.getNumPoints()+1;
//			}
//		}
//	}
//	if(i==ln2.getNumPoints())
//	{
//		ci.covered = false;
//	}
//	i = ln2.getNumPoints()-1;
//	for( ; i>0 && ci.covered; i--)
//	{
//		geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(i), ln2.getCoordinateN(i-1));
//		geos::geom::Coordinate tmp = ln2.getCoordinateN(i);
//		double angle_sg_ln2 = sg_ln2.angle();
//		int j = ln1.getNumPoints() - 1;
//		for( ; j>0; j--)
//		{
//			geos::geom::Coordinate prjc = geos::geom::Coordinate();
//			geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j), ln1.getCoordinateN(j-1));
//			double angle_sg_ln1 = sg_ln1.angle();
//			sg_ln1.project(tmp, prjc);
//			double d_angle1 = fabs(angle_sg_ln1 - angle_sg_ln2);
//			double d_angle2 = 0;
//			if (i < ln2.getNumPoints()-1)
//			{
//				geos::geom::LineSegment sg_ln2_bef = geos::geom::LineSegment(ln2.getCoordinateN(i+1), ln2.getCoordinateN(i));
//				d_angle2 = fabs(sg_ln2_bef.angle() - angle_sg_ln1);
//			}
//			if (sg_ln1.distance(prjc) < 1e-7 && tmp.distance(prjc) < dis_Threshold && (d_angle1 < angle_Threshold || d_angle2 < angle_Threshold))
//			{
//				if(i == ln1.getNumPoints()-1 || prjc.equals(sg_ln2.p0))
//				{
//					ci.endp_ln2 = tmp;
//					ci.endp_ln2_pos = i;
//					ci.endp_ln2_mappingp = prjc;
//					ci.endp_ln2_mappingp_pos = j-1;
//				}
//				else
//				{
//					geos::geom::LineSegment sg_ln2_tmp = geos::geom::LineSegment(ln2.getCoordinateN(i+1), ln2.getCoordinateN(i));
//					geos::geom::Coordinate endp = Mapping::calMappingPoint(sg_ln2_tmp,sg_ln1);
//					if(sg_ln2_tmp.distance(endp) < 1e-7 && sg_ln2.p0.distance(endp) < dis_Threshold)
//					{
//						ci.endp_ln2 = endp;
//						ci.endp_ln2_pos = i;
//						ci.endp_ln2_mappingp = sg_ln1.p0;
//						ci.endp_ln2_mappingp_pos = j;
//					}
//					else
//					{
//						ci.endp_ln2 = tmp;
//						ci.endp_ln2_pos = i;
//						ci.endp_ln2_mappingp = prjc;
//						ci.endp_ln2_mappingp_pos = j-1;
//					}
//				}
//				j = -1;
//				i = -1;
//			}
//		}
//	}
//	if(i==0)
//	{
//		ci.covered = false;
//	}
//
//	return ci;
}

// ln1与ln2之间点豪斯多夫距离 包括ln1到ln2和ln2到ln1的两个豪斯多夫距离的最大和最小值
WayHaustorffDis Mapping::calHausdorffDis(geos::geom::LineString& ln1, geos::geom::LineString& ln2)
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
	WayHaustorffDis wayHausdorff = WayHaustorffDis();
	wayHausdorff.dis_hausdorff_min = std::min(dis_hausdorff_ln1, dis_hausdorff_ln2);
	wayHausdorff.dis_hausdorff_max = std::max(dis_hausdorff_ln1, dis_hausdorff_ln2);
	return wayHausdorff;
}

// ln1与ln2匹配点间的距离 包括最大／最小／平均
WayMappingDis Mapping::calMappingDis(geos::geom::LineString& ln1, geos::geom::LineString& ln2, double dis_Threshold)
{
	CoverInfo ci = Mapping::calCoverOfLine(ln1, ln2, dis_Threshold);
	int num = 0;
	double dis_avg = 0;
	double dis_min = DBL_MAX;
	double dis_max = 0;
	int j = ci.startp_ln1_mappingp_pos;
	if(0 == j)
		j += 1;
	if(!(ci.startp_ln1.equals(ln1.getCoordinateN(ci.startp_ln1_pos))))
	{
		double d = ci.startp_ln1.distance(ci.startp_ln1_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}
	for(int i=ci.startp_ln1_pos; i<=ci.endp_ln1_pos && i<ln1.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		for( ; j<=ci.endp_ln1_mappingp_pos; j++)
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
	if(!(ci.endp_ln1.equals(ln1.getCoordinateN(ci.endp_ln1_pos))))
	{
		double d = ci.endp_ln1.distance(ci.endp_ln1_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}

	j = ci.startp_ln2_mappingp_pos;
	if(0 == j)
		j += 1;
	if(!(ci.startp_ln2.equals(ln2.getCoordinateN(ci.startp_ln2_pos))))
	{
		double d = ci.startp_ln1.distance(ci.startp_ln1_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}
	for(int i=ci.startp_ln2_pos; i<=ci.endp_ln2_pos && i<ln2.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln2 = ln2.getCoordinateN(i);
		for( ; j<=ci.endp_ln2_mappingp_pos && j<ln1.getNumPoints(); j++)
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
	if(!(ci.endp_ln2.equals(ln2.getCoordinateN(ci.endp_ln2_pos))))
	{
		double d = ci.endp_ln2.distance(ci.endp_ln2_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}

	WayMappingDis wayMappingDis = WayMappingDis();
	wayMappingDis.dis_avg = dis_avg / num;
	wayMappingDis.dis_min = dis_min;
	wayMappingDis.dis_max = dis_max;
	return wayMappingDis;
}

// ln1和ln2之间的豪斯多夫距离和匹配距离
WayFeature Mapping::calMappingFeature(geos::geom::LineString& ln1, geos::geom::LineString& ln2, double dis_Threshold)
{
	WayFeature wayFeature = WayFeature();
	CoverInfo ci = Mapping::calCoverOfLine(ln1, ln2, dis_Threshold);
	int num = 0;
	double dis_avg = 0;
	double dis_min = DBL_MAX;
	double dis_max = 0;
	double dis_hausdorff_ln1 = 0;
	double dis_hausdorff_ln2 = 0;
	// 计算ln1起点到匹配起点间到Hausdorff距离
	for(int i=0; i<ci.startp_ln1_pos; i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		if(ci.startp_ln1_mappingp_pos == 0)
		{
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, pt_ln1.distance(ln2.getCoordinateN(0)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=1; j<=ci.startp_ln1_mappingp_pos; j++)
			{
				geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
				d = std::min(d, sg_ln2.distance(pt_ln1));
			}
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, d);
		}
	}
	// 计算ln1第一个匹配点到最后一个匹配点间的Hausdorff距离和匹配距离
	int j = ci.startp_ln1_mappingp_pos;
	if(0 == j)
		j += 1;
	if(!(ci.startp_ln1.equals(ln1.getCoordinateN(ci.startp_ln1_pos))))
	{
		double d = ci.startp_ln1.distance(ci.startp_ln1_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);

	}
	for(int i=ci.startp_ln1_pos; i<=ci.endp_ln1_pos && i<ln1.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		for( ; j<=ci.endp_ln1_mappingp_pos && j<ln2.getNumPoints(); j++)
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
	if(!(ci.endp_ln1.equals(ln1.getCoordinateN(ci.endp_ln1_pos))))
	{
		double d = ci.endp_ln1.distance(ci.endp_ln1_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}
	// 计算最后一个匹配点到最后一个点间到Hausdorff距离
	for(int i=ci.endp_ln1_pos+1; i<ln1.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln1 = ln1.getCoordinateN(i);
		if(ci.endp_ln1_mappingp_pos == ln2.getNumPoints()-1)
		{
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, pt_ln1.distance(ln2.getCoordinateN(ln2.getNumPoints()-1)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=ci.endp_ln1_mappingp_pos+1; j<ln2.getNumPoints(); j++)
			{
				geos::geom::LineSegment sg_ln2 = geos::geom::LineSegment(ln2.getCoordinateN(j-1), ln2.getCoordinateN(j));
				d = std::min(d, sg_ln2.distance(pt_ln1));
			}
			dis_hausdorff_ln1 = std::max(dis_hausdorff_ln1, d);
		}
	}

	// 同上，ln2到ln1
	j = ci.startp_ln2_mappingp_pos;
	if(0 == j)
		j += 1;
	for(int i=0; i<ci.startp_ln2_pos; i++)
	{
		geos::geom::Coordinate pt_ln2 = ln1.getCoordinateN(i);
		if(ci.startp_ln2_mappingp_pos == 0)
		{
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, pt_ln2.distance(ln1.getCoordinateN(0)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=1; j<=ci.startp_ln2_mappingp_pos; j++)
			{
				geos::geom::LineSegment sg_ln1 = geos::geom::LineSegment(ln1.getCoordinateN(j-1), ln1.getCoordinateN(j));
				d = std::min(d, sg_ln1.distance(pt_ln2));
			}
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, d);
		}
	}

	if(!(ci.startp_ln2.equals(ln2.getCoordinateN(ci.startp_ln2_pos))))
	{
		double d = ci.startp_ln1.distance(ci.startp_ln1_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
		std::cout << "--" + std::to_string(d) << std::endl;
	}
	for(int i=ci.startp_ln2_pos; i<=ci.endp_ln2_pos && i<ln2.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln2 = ln2.getCoordinateN(i);
		for( ; j<=ci.endp_ln2_mappingp_pos && j<ln1.getNumPoints(); j++)
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
	if(!(ci.endp_ln2.equals(ln2.getCoordinateN(ci.endp_ln2_pos))))
	{
		double d = ci.endp_ln2.distance(ci.endp_ln2_mappingp);
		dis_avg += d;
		num += 1;
		dis_min = std::min(dis_min, d);
		dis_max = std::max(dis_max, d);
	}

	for(int i=ci.endp_ln2_pos+1; i<ln2.getNumPoints(); i++)
	{
		geos::geom::Coordinate pt_ln2 = ln2.getCoordinateN(i);
		if(ci.endp_ln2_mappingp_pos == ln1.getNumPoints()-1)
		{
			dis_hausdorff_ln2 = std::max(dis_hausdorff_ln2, pt_ln2.distance(ln1.getCoordinateN(ln1.getNumPoints()-1)));
		}
		else
		{
			double d = DBL_MAX;
			for(int j=ci.endp_ln2_mappingp_pos+1; j<ln1.getNumPoints(); j++)
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
	wayFeature.dis_hausdorff_min = std::min(dis_hausdorff_ln1, dis_hausdorff_ln2);
	wayFeature.dis_hausdorff_max = std::max(dis_hausdorff_ln1, dis_hausdorff_ln2);
	wayFeature.dis_hausdorff_avg = (dis_hausdorff_ln1 + dis_hausdorff_ln2) / 2;
	return wayFeature;
}
