/*
 * feature_probability_calculator.cpp
 *
 *  Created on: 2018年5月10日
 *      Author: yinweijun
 */

#include "matching_probability_calculator.h"

#include <cmath>

const double CHANGED = 300;

MatchingProbabilityCalculator::MatchingProbabilityCalculator(
		WayOperator* way_operator1, WayOperator* way_operator2,
		const MapWayMatchingPair& map_way_matching_pair)
:  way_operator1_(way_operator1)
, way_operator2_(way_operator2)
, map_way_matching_pair_(map_way_matching_pair)
{
	init();
}

MatchingProbabilityCalculator::~MatchingProbabilityCalculator()
{
}

void MatchingProbabilityCalculator::calculator()
{
	update();
}

void MatchingProbabilityCalculator::init()
{
	int rows = way_operator1_->map_way_osm_id_.size();
	int cols = way_operator2_->map_way_osm_id_.size();
	matrix_confidence_.resize(rows, cols);
	matrix_support_.resize(rows, cols);
	matrix_hausdorff_dis0_.resize(rows, cols);
	matrix_hausdorff_dis1_.resize(rows, cols);
	matrix_hausdorff_dir0_.resize(rows, cols);
	matrix_hausdorff_dir1_.resize(rows, cols);
	MatchingProbabilityCalculator::initEigen(ids_main_, ids_mapped_);
}

void MatchingProbabilityCalculator::initEigen(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
{
	int num_rows(0);
	int num_cols(0);
	MapWayMatchingPair::iterator it_wmp = map_way_matching_pair_.begin();
	for (; it_wmp != map_way_matching_pair_.end(); ++it_wmp)
	{
		PWayMatchingPair pWayMatchingPair = it_wmp->second;
		PWay way_main = pWayMatchingPair->main_way;
		PWay way_mapped = pWayMatchingPair->mapped_way;

		int index_main(num_rows);
		int index_mapped(num_cols);
		if(ids_main.find(way_main->osm_id)!=ids_main.end())
		{
			index_main = ids_main[way_main->osm_id];
		}
		else
		{
			ids_main.insert(make_pair(way_main->osm_id, num_rows));
			num_rows += 1;
		}
		if(ids_mapped.find(way_mapped->osm_id)!=ids_mapped.end())
		{
			index_mapped = ids_mapped[way_mapped->osm_id];
		}
		else
		{
			ids_mapped.insert(make_pair(way_mapped->osm_id, num_cols));
			num_cols += 1;
		}
//		std::cout << "+++++++++++++++++++++++++++++++++++++++++" << std::endl;
//		std::cout << "	index_main:	" << index_main << "	index_mapped:	" << index_mapped << std::endl;
//		std::cout << "	SIZE:	row:" << matrix_confidence_.rows() << ", cols" << matrix_confidence_.cols() << std::endl;
//		std::cout << "		confidence: " << pWayMatchingPair->confidence << ",  dis_hausdorff_ln1: " <<  pWayMatchingPair->dis_hausdorff_ln1 << " ,  angle_hausdorff_ln1: " << pWayMatchingPair->angle_hausdorff_ln1 << std::endl;
		matrix_confidence_(index_main, index_mapped) = pWayMatchingPair->confidence;
		matrix_hausdorff_dis0_(index_main, index_mapped) = pWayMatchingPair->dis_hausdorff_ln1;
		matrix_hausdorff_dis1_(index_main, index_mapped) = pWayMatchingPair->dis_hausdorff_ln2;
		matrix_hausdorff_dir0_(index_main, index_mapped) = pWayMatchingPair->angle_hausdorff_ln1;
		matrix_hausdorff_dir1_(index_main, index_mapped) = pWayMatchingPair->angle_hausdorff_ln2;

	}
}

void MatchingProbabilityCalculator::calQIJ(double& q, PWayMatchingPair& pWayMatchingPair, std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped, bool last)
{
	bool debug(false);
	PWay way_main = pWayMatchingPair->main_way;
	PWay way_mapped = pWayMatchingPair->mapped_way;
	int index_way_main = ids_main[way_main->osm_id];
	int index_way_mapped = ids_mapped[way_mapped->osm_id];

	double dis_hausdorff0(0.);
	double dis_hausdorff1(0.);
	std::string key(way_main->osm_id + "$" + way_mapped->osm_id);
	MapWayMatchingPair::iterator it_find = map_way_matching_pair_.find(key);
	if (it_find != map_way_matching_pair_.end())
	{
		dis_hausdorff0 = it_find->second->dis_mapping_sta;
		dis_hausdorff1 = it_find->second->dis_mapping_end;
	}

	double angle_start_main = way_main->sta_angle;
	double angle_end_main = way_main->end_angle;
	double angle_start_mapped = way_mapped->sta_angle;
	double angle_end_mapped = way_mapped->end_angle;

	std::vector<PWay> main;
	std::vector<PWay> mapped;

	if(last)
	{
		if (pWayMatchingPair->same_direction)
		{
			main = way_main->vec_last_way;
			mapped = way_mapped->vec_last_way;
		}
		else
		{
			main = way_main->vec_last_way;
			mapped = way_mapped->vec_next_way;
		}
	}
	else
	{
		if (pWayMatchingPair->same_direction)
		{
			main = way_main->vec_next_way;
			mapped = way_mapped->vec_next_way;
		}
		else
		{
			main = way_main->vec_next_way;
			mapped = way_mapped->vec_last_way;
		}
	}

	if (way_main->osm_id == "w1712502152726" && way_mapped->osm_id == "w1709001738553")
	{
		debug = true;
	}
	if(debug)
	{
		std::cout << "==========================IJ=========================================" << std::endl;
		std::cout << "	Confidence:  " << pWayMatchingPair->confidence << ",	dir1:  " << pWayMatchingPair->angle_hausdorff_ln1 << ",	dir2:  " << pWayMatchingPair->angle_hausdorff_ln2 << ",	dis1:  " << dis_hausdorff0 << ",	dis2:  " << dis_hausdorff1 << std::endl;
		std::cout << "	main_angle_sta:  " << angle_start_main << ",	main_angle_end: " << angle_end_main << std::endl;
		std::cout << "	mapped_angle_sta:  " << angle_start_mapped << ",	mapped_angle_end: " << angle_end_mapped << std::endl;
	}

	double q1(0.);
	double q2(0.);

	// 没有last
	if(main.size() < 1 || mapped.size() < 1)
	{
		if(main.size() < 1 && mapped.size() < 1)
		{
			q1 = 1.0;
			q2 = 1.0;
		}
		else if(main.size() < 1)
		{
			std::vector<PWay>::iterator it_mapped = mapped.begin();
			double p_last_mapped(0.);
			int num_p_last_mapped(0);
			for(; it_mapped != mapped.end(); ++it_mapped)
			{
				PWay way_last_mapped = *it_mapped;
				std::map<std::string, int>::iterator it_find_last_mapped = ids_mapped.find(way_last_mapped->osm_id);
				if (it_find_last_mapped == ids_mapped.end())
					continue;
				int index_way_last_mapped = ids_mapped[way_last_mapped->osm_id];
				if(matrix_confidence_.col(index_way_last_mapped).maxCoeff() < 0.1)
					continue;
				if(matrix_confidence_(index_way_main, index_way_last_mapped) > 1e-7)
				{
					p_last_mapped += matrix_confidence_(index_way_main, index_way_last_mapped);
					num_p_last_mapped += 1;
				}
				else
				{
					p_last_mapped += 1.0 - matrix_confidence_.col(index_way_last_mapped).maxCoeff();
					num_p_last_mapped += 1;
				}
			}
			if (num_p_last_mapped > 1)
				p_last_mapped = p_last_mapped / num_p_last_mapped;
			q1 = p_last_mapped;
			q2 = p_last_mapped;
		}
		else
		{
			std::vector<PWay>::iterator it_main = main.begin();
			double p_last_main(0.);
			int num_p_last_main(0);
			for(; it_main != main.end(); ++it_main)
			{
				PWay way_last_main = *it_main;
				std::map<std::string, int>::iterator it_find_last_main = ids_main.find(way_last_main->osm_id);
				if (it_find_last_main == ids_main.end())
					continue;
				int index_way_last_main = ids_main[way_last_main->osm_id];
				if(matrix_confidence_.row(index_way_last_main).maxCoeff() < 0.1)
					continue;
				if (matrix_confidence_(index_way_last_main, index_way_mapped) > 1e-7)
				{
					p_last_main += matrix_confidence_(index_way_last_main, index_way_mapped);
					num_p_last_main += 1;
				}
				else
				{
					p_last_main += 1.0 - matrix_confidence_.row(index_way_last_main).maxCoeff();
					num_p_last_main += 1;
				}
			}
			if (num_p_last_main > 1)
				p_last_main = p_last_main / num_p_last_main;
			q1 = p_last_main;
			q2 = p_last_main;
		}

		q = std::max(q1, q2);
		if(debug)
		{
			std::cout << "\n	IJ Result:  " << q << std::endl;
		}
		return;
	}

	// q1
	int num_q1(0);
	std::vector<PWay>::iterator it_main = main.begin();
	int num_add(0);
	for(; it_main != main.end(); ++it_main)
	{
		PWay way_last_main = *it_main;
		// TODO 1:0
		std::map<std::string, int>::iterator it_find_last_main = ids_main.find(way_last_main->osm_id);
		if (it_find_last_main == ids_main.end())
		{
			num_add += 1;
			continue;
		}
		int index_way_last_main = ids_main[way_last_main->osm_id];

		std::vector<PWay>::iterator it_last_mapped = mapped.begin();
		double cp(0.);
		for(; it_last_mapped != mapped.end(); ++it_last_mapped)
		{
			PWay way_last_mapped = *it_last_mapped;
			std::map<std::string, int>::iterator it_find_last_mapped = ids_mapped.find(way_last_mapped->osm_id);
			if (it_find_last_mapped == ids_mapped.end())
				continue;

			int index_way_last_mapped = ids_mapped[way_last_mapped->osm_id];

			// i,j,a,b
			double confidence = matrix_confidence_(index_way_last_main, index_way_last_mapped);
			if(confidence > 1e-7 )
			{
				if(debug)
					std::cout << "	Confidence: " << confidence << ",	maxCoeff: " <<  matrix_confidence_.row(index_way_last_main).maxCoeff() << std::endl;
				double d1(0.);
				double d2(0.);
				double rate(0.);
				std::string key(way_last_main->osm_id + "$" + way_last_mapped->osm_id);
				MapWayMatchingPair::iterator it_find = map_way_matching_pair_.find(key);
				if (it_find != map_way_matching_pair_.end())
				{
					d1 = it_find->second->dis_mapping_sta;
					d2 = it_find->second->dis_mapping_end;
					if(it_find->second->length_mapping_ln2 < 1e-7)
						continue;
					rate = fabs(1.0 - fabs(it_find->second->length_mapping_ln1 - it_find->second->length_mapping_ln2) / std::max(it_find->second->length_mapping_ln1, it_find->second->length_mapping_ln2));
				}
				else
					continue;
				double q_dis = 0.0;
				double q_dir = 0.0;
				if((dis_hausdorff0 + dis_hausdorff1) > 1e-7)
				{
					q_dis = fabs((d1 + d2) / 2 - (dis_hausdorff0 + dis_hausdorff1) / 2);
				}
				double d_dir(0.);
				if(last)
				{
					if (pWayMatchingPair->same_direction)
						d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->end_angle, angle_start_main), GeosTool::calculate_include_angle(way_last_mapped->end_angle, angle_start_mapped)));
					else
						d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->end_angle, angle_start_main), GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped)));

					if(debug)
					{
						std::cout << "	last_main_end:  " << way_last_main->end_angle << "	main_start:  " << angle_start_main << "	dir1:  " << GeosTool::calculate_include_angle(way_last_main->end_angle, angle_start_main) << std::endl;
						std::cout << "	last_mapped_sta:  " << way_last_mapped->sta_angle << "	mapped_end:  " << angle_end_mapped << "	dir1:  " << GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped) << std::endl;
					}
				}
				else
				{
					if (pWayMatchingPair->same_direction)
						d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->sta_angle, angle_end_main), GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped)));
					else
						d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->sta_angle, angle_end_main), GeosTool::calculate_include_angle(way_last_mapped->end_angle, angle_start_mapped)));
				}
				q_dir = d_dir;

				double a = -0.2*std::log(2);
				double b = 4 * std::log(2);
				double dis = 1.0 / (1.0 + std::pow(M_E, -(a * q_dis + b)));
				a = -8*std::log(2);
				b = 6 * std::log(2);
				double dir = 1.0 / (1.0 + std::pow(M_E, -(a * q_dir + b)));

				if (q_dis < 1e-7)
					dis = 1.0;
				if (q_dir < 1e-7)
					dir = 1.0;

				double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * rate;

				if(debug)
				{
					std::cout << "-----------------i,j,a,b-----------------------------" << std::endl;
					std::cout << "		main_last_way: " << way_last_main->osm_id << ",	mapped_last_way: " << way_last_mapped->osm_id << std::endl;
					std::cout << "		main_last_way_sta_angle: " << way_last_main->sta_angle << ",	main_last_way_end_angle: " << way_last_main->end_angle << std::endl;
					std::cout << "		mapped_last_way_sta_angle: " << way_last_mapped->sta_angle << ",	mapped_last_way_end_angle: " << way_last_mapped->end_angle << std::endl;
					std::cout << "		dis1: " << d1 << ",	dis2: " << d2 << std::endl;
					std::cout << "		d_dir: " << d_dir << std::endl;
					std::cout << "		q_dis: " << q_dis << ",	q_dir: " << q_dir << std::endl;
					std::cout << "		length_mapping_ln2: " << it_find->second->length_mapping_ln2 << ",	length_mapping_ln1: " << it_find->second->length_mapping_ln1 << std::endl;
					std::cout << "		rate:  " << rate  << std::endl;
					std::cout << "		c: " << c << ",	cp(q1): " << c * confidence << std::endl;
				}

				cp = std::max(cp, c * confidence);
			}
		}
		// i,j,a.j
		double confidence_j = matrix_confidence_(index_way_last_main, index_way_mapped);
		if(confidence_j > 1e-7 )
		{
			if(debug)
				std::cout << "	Confidence: " << confidence_j << ",	maxCoeff: " <<  matrix_confidence_.row(index_way_last_main).maxCoeff() << std::endl;
			double d1(0.);
			double d2(0.);
			double rate(0.);
			std::string key(way_last_main->osm_id + "$" + way_mapped->osm_id);
			MapWayMatchingPair::iterator it_find = map_way_matching_pair_.find(key);
			if (it_find != map_way_matching_pair_.end())
			{
				d1 = it_find->second->dis_mapping_sta;
				d2 = it_find->second->dis_mapping_end;
				if(it_find->second->length_mapping_ln2 < 1e-7)
				{
					if(cp > 1e-7)
					{
						if (matrix_confidence_.row(index_way_last_main).maxCoeff() > 0.1)
						{
							num_q1 += 1;
							q1 += cp;
						}
						else
						{
							num_add += 1;
						}
					}
					else
						num_add += 1;
					continue;
				}
				rate = 1.0 - fabs(it_find->second->length_mapping_ln1 - it_find->second->length_mapping_ln2) / std::max(it_find->second->length_mapping_ln1, it_find->second->length_mapping_ln2);
			}
			else
			{
				if(cp > 1e-7)
				{
					if(matrix_confidence_.row(index_way_last_main).maxCoeff() > 0.1)
					{
						num_q1 += 1;
						q1 += cp;
					}
					else
						num_add += 1;
				}
				else
					num_add += 1;
				continue;
			}
			double q_dis = 0.0;
			double q_dir = 0.0;
			if((dis_hausdorff0 + dis_hausdorff1) > 1e-7)
			{
				q_dis = fabs((d1 + d2) / 2 - (dis_hausdorff0 + dis_hausdorff1) / 2);
			}
			double d_dir(0.);
//			if(last)
//			{
//				d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->end_angle, angle_start_main), 0));
//			}
//			else
//			{
//				d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->sta_angle, angle_end_main), 0));
//			}
			q_dir = d_dir;

			double a = -0.2*std::log(2);
			double b = 4 * std::log(2);
			double dis = 1.0 / (1.0 + std::pow(M_E, -(a * q_dis + b)));
			a = -8*std::log(2);
			b = 6 * std::log(2);
			double dir = 1.0 / (1.0 + std::pow(M_E, -(a * q_dir + b)));

			if (q_dis < 1e-7)
				dis = 1.0;
			if (q_dir < 1e-7)
				dir = 1.0;

			double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * rate;

			if(debug)
			{
				std::cout << "---------------------i,j,a,j-------------------------" << std::endl;
				std::cout << "		main_last_way: " << way_last_main->osm_id << ",	mapped_last_way: " << way_mapped->osm_id << std::endl;
				std::cout << "		main_last_way_sta_angle: " << way_last_main->sta_angle << ",	main_last_way_end_angle: " << way_last_main->end_angle << std::endl;
				std::cout << "		mapped_last_way_sta_angle: " << way_mapped->sta_angle << ",	mapped_last_way_end_angle: " << way_mapped->end_angle << std::endl;
				std::cout << "		dis1: " << d1 << ",	dis2: " << d2 << std::endl;
				std::cout << "		d_dir: " << d_dir << std::endl;
				std::cout << "		q_dis: " << q_dis << ",	q_dir: " << q_dir << std::endl;
				std::cout << "		length_mapping_ln2: " << it_find->second->length_mapping_ln2 << ",	length_mapping_ln1: " << it_find->second->length_mapping_ln1 << std::endl;
				std::cout << "		rate:  " << rate  << std::endl;
				std::cout << "		c: " << c << ",	cp(q1): " << c * confidence_j << std::endl;
			}

			cp = std::max(cp, c * confidence_j);
		}
		if(cp > 1e-7)
		{
			if(matrix_confidence_.row(index_way_last_main).maxCoeff() > 0.1)
			{
				num_q1 += 1;
				q1 += cp;
			}
			else
				num_add += 1;
		}
		else
		{
			if(matrix_confidence_.row(index_way_last_main).maxCoeff() < 0.1)
				num_add += 1;
		}
	}
	if(main.size() - num_add > 1e-7)
		q1 = q1 / (main.size() - num_add);
	else
		q1 = matrix_confidence_(index_way_main, index_way_mapped);

	//q2
	// i,j,i,b
	std::vector<PWay>::iterator it_last_mapped = mapped.begin();
	for(; it_last_mapped != mapped.end(); ++it_last_mapped)
	{
		PWay way_last_mapped = *it_last_mapped;
		std::map<std::string, int>::iterator it_find_last_mapped = ids_mapped.find(way_last_mapped->osm_id);
		if (it_find_last_mapped == ids_mapped.end())
			continue;
		int index_way_last_mapped = ids_mapped[way_last_mapped->osm_id];
		double confidence_i = matrix_confidence_(index_way_main, index_way_last_mapped);
		if(confidence_i > 1e-7 )
		{
			double d1(0.);
			double d2(0.);
			double rate(0.);
			std::string key(way_main->osm_id + "$" + way_last_mapped->osm_id);
			MapWayMatchingPair::iterator it_find = map_way_matching_pair_.find(key);
			if (it_find != map_way_matching_pair_.end())
			{
				d1 = it_find->second->dis_mapping_sta;
				d2 = it_find->second->dis_mapping_end;
				if(it_find->second->length_mapping_ln2 < 1e-7)
					continue;
				rate = fabs(1.0 - fabs(it_find->second->length_mapping_ln1 - it_find->second->length_mapping_ln2) / std::max(it_find->second->length_mapping_ln1, it_find->second->length_mapping_ln2));
			}
			else
				continue;

			double q_dis = 0.0;
			double q_dir = 0.0;
			if((dis_hausdorff0 + dis_hausdorff1) > 1e-7)
			{
				q_dis = fabs((d1 + d2) / 2 - (dis_hausdorff0 + dis_hausdorff1) / 2);
			}
			double d_dir(0.);
//			if(last)
//			{
//				if (pWayMatchingPair->same_direction)
//					d_dir = fabs(GeosTool::calculate_include_angle(0, GeosTool::calculate_include_angle(way_last_mapped->end_angle, angle_start_mapped)));
//				else
//					d_dir = fabs(GeosTool::calculate_include_angle(0, GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped)));
//			}
//			else
//			{
//				if (pWayMatchingPair->same_direction)
//					d_dir = fabs(GeosTool::calculate_include_angle(0, GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped)));
//				else
//					d_dir = fabs(GeosTool::calculate_include_angle(0, GeosTool::calculate_include_angle(way_last_mapped->end_angle, angle_start_mapped)));
//			}
			q_dir = d_dir;

			double a = -0.2*std::log(2);
			double b = 4 * std::log(2);
			double dis = 1.0 / (1.0 + std::pow(M_E, -(a * q_dis + b)));
			a = -8*std::log(2);
			b = 6 * std::log(2);
			double dir = 1.0 / (1.0 + std::pow(M_E, -(a * q_dir + b)));

			if (q_dis < 1e-7)
				dis = 1.0;
			if (q_dir < 1e-7)
				dir = 1.0;

			double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * rate;

			if(debug)
			{
				std::cout << "---------------------i,j,i,b-------------------------" << std::endl;
				std::cout << "		main_last_way: " << way_main->osm_id << ",	mapped_last_way: " << way_last_mapped->osm_id << std::endl;
				std::cout << "		main_last_way_sta_angle: " << way_main->sta_angle << ",	main_last_way_end_angle: " << way_main->end_angle << std::endl;
				std::cout << "		mapped_last_way_sta_angle: " << way_last_mapped->sta_angle << ",	mapped_last_way_end_angle: " << way_last_mapped->end_angle << std::endl;
				std::cout << "		dis1: " << d1 << ",	dis2: " << d2 << std::endl;
				std::cout << "		d_dir: " << d_dir << std::endl;
				std::cout << "		q_dis: " << q_dis << ",	q_dir: " << q_dir << std::endl;
				std::cout << "		length_mapping_ln2: " << it_find->second->length_mapping_ln2 << ",	length_mapping_ln1: " << it_find->second->length_mapping_ln1 << std::endl;
				std::cout << "		rate:  " << rate  << std::endl;
				std::cout << "		c: " << c << ",	cp(q2): " << c * confidence_i << std::endl;
			}

			q2 = std::max(q2, c * confidence_i);
		}
	}

//	if(q1 > 1e-7 || q2 > 1e-7)
//	{
//		if(q1 > 1e-7 && q2 > 1e-7)
//		{
//			q = (q1 + q2) / (main.size() - num_add + 1);
//		}
//		else if(q1 > 1e-7)
//		{
//			if(main.size() - num_add > 1e-7)
//			{
//				q = q1 / (main.size() - num_add);
//			}
//			else
//			{
//				q = matrix_confidence_(index_way_main, index_way_mapped);
//			}
//		}
//		else
//		{
//			q = q2;
//		}
//	}
//	else
//		q = 0.0;
	q = std::max(q1, q2);
	if(debug)
	{
		std::cout << "JI: num_q1: " << num_q1 << std::endl;
		std::cout << "q1: " << q1 << ", q2: " << q2 << std::endl;
		std::cout << "\n	IJ Result:  " << q << ",	num_added: " << num_add << std::endl;
	}
}

void MatchingProbabilityCalculator::calQJI(double& q, PWayMatchingPair& pWayMatchingPair, std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped, bool last)
{
	bool debug(false);
	PWay way_main = pWayMatchingPair->main_way;
	PWay way_mapped = pWayMatchingPair->mapped_way;
	int index_way_main = ids_main[way_main->osm_id];
	int index_way_mapped = ids_mapped[way_mapped->osm_id];

	double dis_hausdorff0(0.);
	double dis_hausdorff1(0.);
	std::string key(way_main->osm_id + "$" + way_mapped->osm_id);
	MapWayMatchingPair::iterator it_find = map_way_matching_pair_.find(key);
	if (it_find != map_way_matching_pair_.end())
	{
		dis_hausdorff0 = it_find->second->dis_mapping_sta;
		dis_hausdorff1 = it_find->second->dis_mapping_end;
	}

	double angle_start_main = way_main->sta_angle;
	double angle_end_main = way_main->end_angle;
	double angle_start_mapped = way_mapped->sta_angle;
	double angle_end_mapped = way_mapped->end_angle;

	std::vector<PWay> main;
	std::vector<PWay> mapped;

	if(last)
	{
		if (pWayMatchingPair->same_direction)
		{
			main = way_main->vec_last_way;
			mapped = way_mapped->vec_last_way;
		}
		else
		{
			main = way_main->vec_last_way;
			mapped = way_mapped->vec_next_way;
		}
	}
	else
	{
		if (pWayMatchingPair->same_direction)
		{
			main = way_main->vec_next_way;
			mapped = way_mapped->vec_next_way;
		}
		else
		{
			main = way_main->vec_next_way;
			mapped = way_mapped->vec_last_way;
		}
	}
	if (way_main->osm_id == "w1712502152726" && way_mapped->osm_id == "w1709001738553")
	{
		debug = true;
	}
	if(debug)
	{
		std::cout << "==========================JI=========================================" << std::endl;
		std::cout << "	Confidence:  " << pWayMatchingPair->confidence << ",	dir1:  " << pWayMatchingPair->angle_hausdorff_ln1 << ",	dir2:  " << pWayMatchingPair->angle_hausdorff_ln2 << ",	dis1:  " << dis_hausdorff0 << ",	dis2:  " << dis_hausdorff1 << std::endl;
		std::cout << "	main_angle_sta:  " << angle_start_main << ",	main_angle_end: " << angle_end_main << std::endl;
		std::cout << "	mapped_angle_sta:  " << angle_start_mapped << ",	mapped_angle_end: " << angle_end_mapped << std::endl;
	}

	double q1(0.);
	double q2(0.);

	// 没有last
	if(main.size() < 1 || mapped.size() < 1)
	{
		if(main.size() < 1 && mapped.size() < 1)
		{
			q1 = 1.0;
			q2 = 1.0;
		}
		else if(main.size() < 1)
		{
			std::vector<PWay>::iterator it_mapped = mapped.begin();
			double p_last_mapped(0.);
			int num_p_last_mapped(0);
			for(; it_mapped != mapped.end(); ++it_mapped)
			{
				PWay way_last_mapped = *it_mapped;
				std::map<std::string, int>::iterator it_find_last_mapped = ids_mapped.find(way_last_mapped->osm_id);
				if (it_find_last_mapped == ids_mapped.end())
					continue;
				int index_way_last_mapped = ids_mapped[way_last_mapped->osm_id];
				if(matrix_confidence_.col(index_way_last_mapped).maxCoeff() < 0.1)
					continue;
				if (matrix_confidence_(index_way_main, index_way_last_mapped) > 1e-7)
				{
					p_last_mapped += matrix_confidence_(index_way_main, index_way_last_mapped);
					num_p_last_mapped += 1;
				}
				else
				{
					p_last_mapped += 1.0 - matrix_confidence_.col(index_way_last_mapped).maxCoeff();
					num_p_last_mapped += 1;
				}
			}
			if (num_p_last_mapped > 1)
				p_last_mapped = p_last_mapped / num_p_last_mapped;
			q1 = p_last_mapped;
			q2 = p_last_mapped;
		}
		else
		{
			std::vector<PWay>::iterator it_main = main.begin();
			double p_last_main(0.);
			int num_p_last_main(0);
			for(; it_main != main.end(); ++it_main)
			{
				PWay way_last_main = *it_main;
				std::map<std::string, int>::iterator it_find_last_main = ids_main.find(way_last_main->osm_id);
				if (it_find_last_main == ids_main.end())
					continue;
				int index_way_last_main = ids_main[way_last_main->osm_id];
				if(matrix_confidence_.row(index_way_last_main).maxCoeff() < 0.1)
					continue;
				if (matrix_confidence_(index_way_last_main, index_way_mapped) > 1e-7)
				{
					p_last_main += matrix_confidence_(index_way_last_main, index_way_mapped);
					num_p_last_main += 1;
				}
				else
				{
					p_last_main += 1.0 - matrix_confidence_.row(index_way_last_main).maxCoeff();
					num_p_last_main += 1;
				}
			}
			if (num_p_last_main > 1)
				p_last_main = p_last_main / num_p_last_main;
			q1 = p_last_main;
			q2 = p_last_main;
		}

		q = std::max(q1, q2);
		if(debug)
		{
			std::cout << "\n	JI Result:  " << q << std::endl;
		}
		return;
	}

	// q1
	int num_q1(0);
	std::vector<PWay>::iterator it_mapped = mapped.begin();
	int num_add(0);
	for(; it_mapped != mapped.end(); ++it_mapped)
	{
		PWay way_last_mapped = *it_mapped;
		std::map<std::string, int>::iterator it_find_last_mapped = ids_mapped.begin();
		if (it_find_last_mapped == ids_mapped.end())
		{
			num_add += 1;
			continue;
		}
		int index_way_last_mapped = ids_mapped[way_last_mapped->osm_id];

		std::vector<PWay>::iterator it_last_main = main.begin();
		double cp(0.);
		for(; it_last_main != main.end(); ++it_last_main)
		{
			PWay way_last_main = *it_last_main;
			std::map<std::string, int>::iterator it_find_last_main = ids_main.begin();
			if (it_find_last_main == ids_main.end())
				continue;
			int index_way_last_main = ids_main[way_last_main->osm_id];
			// j,i,b,a
			double confidence = matrix_confidence_(index_way_last_main, index_way_last_mapped);
			if(confidence > 1e-7 )
			{
				if(debug)
					std::cout << "	Confidence: " << confidence << ",	maxCoeff: " <<  matrix_confidence_.col(index_way_last_mapped).maxCoeff() << std::endl;
				double d1(0.);
				double d2(0.);
				double rate(0.);
				std::string key(way_last_main->osm_id + "$" + way_last_mapped->osm_id);
				MapWayMatchingPair::iterator it_find = map_way_matching_pair_.find(key);
				if (it_find != map_way_matching_pair_.end())
				{
					d1 = it_find->second->dis_mapping_sta;
					d2 = it_find->second->dis_mapping_end;
					if(it_find->second->length_mapping_ln1 < 1e-7)
						continue;
					rate = fabs(1.0 - fabs(it_find->second->length_mapping_ln1 - it_find->second->length_mapping_ln2) / std::max(it_find->second->length_mapping_ln1, it_find->second->length_mapping_ln2));
				}
				else
					continue;
				double q_dis = 0.0;
				double q_dir = 0.0;
				if((dis_hausdorff0 + dis_hausdorff1) > 1e-7)
				{
					q_dis = fabs((d1 + d2) / 2 - (dis_hausdorff0 + dis_hausdorff1) / 2);
				}
				double d_dir(0.);
				if(last)
				{
					if (pWayMatchingPair->same_direction)
						d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->end_angle, angle_start_main), GeosTool::calculate_include_angle(way_last_mapped->end_angle, angle_start_mapped)));
					else
						d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->end_angle, angle_start_main), GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped)));
				}
				else
				{
					if (pWayMatchingPair->same_direction)
						d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->sta_angle, angle_end_main), GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped)));
					else
						d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_main->sta_angle, angle_end_main), GeosTool::calculate_include_angle(way_last_mapped->end_angle, angle_start_mapped)));
				}
				q_dir = d_dir;

				double a = -0.2*std::log(2);
				double b = 4 * std::log(2);
				double dis = 1.0 / (1.0 + std::pow(M_E, -(a * q_dis + b)));
				a = -8*std::log(2);
				b = 6 * std::log(2);
				double dir = 1.0 / (1.0 + std::pow(M_E, -(a * q_dir + b)));

				if (q_dis < 1e-7)
					dis = 1.0;
				if (q_dir < 1e-7)
					dir = 1.0;

				double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * rate;

				if(debug)
				{
					std::cout << "-----------------j,i,b,a-----------------------------" << std::endl;
					std::cout << "		main_last_way: " << way_last_main->osm_id << ",	mapped_last_way: " << way_last_mapped->osm_id << std::endl;
					std::cout << "		main_last_way_sta_angle: " << way_last_main->sta_angle << ",	main_last_way_end_angle: " << way_last_main->end_angle << std::endl;
					std::cout << "		mapped_last_way_sta_angle: " << way_last_mapped->sta_angle << ",	mapped_last_way_end_angle: " << way_last_mapped->end_angle << std::endl;
					std::cout << "		dis1: " << d1 << ",	dis2: " << d2 << std::endl;
					std::cout << "		d_dir: " << d_dir << std::endl;
					std::cout << "		q_dis: " << q_dis << ",	q_dir: " << q_dir << std::endl;
					std::cout << "		length_mapping_ln2: " << it_find->second->length_mapping_ln2 << ",	length_mapping_ln1: " << it_find->second->length_mapping_ln1 << std::endl;
					std::cout << "		rate:  " << rate  << std::endl;
					std::cout << "		c: " << c << ",	cp(q1): " << c * confidence << std::endl;
				}

				cp = std::max(cp, c * confidence);
			}
		}
		// j,i,b.i
		double confidence_i = matrix_confidence_(index_way_main, index_way_last_mapped);
		if(confidence_i > 1e-7 )
		{
			if(debug)
				std::cout << "	Confidence: " << confidence_i << ",	maxCoeff: " <<  matrix_confidence_.col(index_way_last_mapped).maxCoeff() << std::endl;
			double d1(0.);
			double d2(0.);
			double rate(0.);
			std::string key(way_main->osm_id + "$" + way_last_mapped->osm_id);
			MapWayMatchingPair::iterator it_find = map_way_matching_pair_.find(key);
			if (it_find != map_way_matching_pair_.end())
			{
				d1 = it_find->second->dis_mapping_sta;
				d2 = it_find->second->dis_mapping_end;
				if(it_find->second->length_mapping_ln1 < 1e-7)
				{
					if (cp > 1e-7)
					{
						if(matrix_confidence_.col(index_way_last_mapped).maxCoeff() > 0.1)
						{
							num_q1 += 1;
							q1 += cp;
						}
						else
							num_add += 1;
					}
					else
						num_add += 1;
					continue;
				}
				rate = 1 - fabs(it_find->second->length_mapping_ln2 - it_find->second->length_mapping_ln1) / std::max(it_find->second->length_mapping_ln1, it_find->second->length_mapping_ln2);
			}
			else
			{
				if (cp > 1e-7)
				{
					if(matrix_confidence_.col(index_way_last_mapped).maxCoeff() > 0.1)
					{
						q1 += cp;
						num_q1 += 1;
					}
				}
				else
					num_add += 1;
				continue;
			}

			double q_dis = 0.0;
			double q_dir = 0.0;
			if((dis_hausdorff0 + dis_hausdorff1) > 1e-7)
			{
				q_dis = fabs((d1 + d2) / 2 - (dis_hausdorff0 + dis_hausdorff1) / 2);
			}
			double d_dir(0.);
//			if(last)
//			{
//				if (pWayMatchingPair->same_direction)
//					d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_mapped->end_angle, angle_start_mapped), 0));
//				else
//					d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped), 0));
//			}
//			else
//			{
//				if (pWayMatchingPair->same_direction)
//					d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_mapped->sta_angle, angle_end_mapped), 0));
//				else
//					d_dir = fabs(GeosTool::calculate_include_angle(GeosTool::calculate_include_angle(way_last_mapped->end_angle, angle_start_mapped), 0));
//			}
			q_dir = d_dir;

			double a = -0.2*std::log(2);
			double b = 4 * std::log(2);
			double dis = 1.0 / (1.0 + std::pow(M_E, -(a * q_dis + b)));
			a = -8*std::log(2);
			b = 6 * std::log(2);
			double dir = 1.0 / (1.0 + std::pow(M_E, -(a * q_dir + b)));

			if (q_dis < 1e-7)
				dis = 1.0;
			if (q_dir < 1e-7)
				dir = 1.0;

			double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * rate;

			if(debug)
			{
				std::cout << "---------------------i,j,a,j-------------------------" << std::endl;
				std::cout << "		main_last_way: " << way_main->osm_id << ",	mapped_last_way: " << way_last_mapped->osm_id << std::endl;
				std::cout << "		main_last_way_sta_angle: " << way_main->sta_angle << ",	main_last_way_end_angle: " << way_main->end_angle << std::endl;
				std::cout << "		mapped_last_way_sta_angle: " << way_last_mapped->sta_angle << ",	mapped_last_way_end_angle: " << way_last_mapped->end_angle << std::endl;
				std::cout << "		dis1: " << d1 << ",	dis2: " << d2 << std::endl;
				std::cout << "		d_dir: " << d_dir << std::endl;
				std::cout << "		q_dis: " << q_dis << ",	q_dir: " << q_dir << std::endl;
				std::cout << "		length_mapping_ln2: " << it_find->second->length_mapping_ln2 << ",	length_mapping_ln1: " << it_find->second->length_mapping_ln1 << std::endl;
				std::cout << "		rate:  " << rate  << std::endl;
				std::cout << "		c: " << c << ",	cp(q2): " << c * confidence_i << std::endl;
			}

			cp = std::max(cp, c * confidence_i);
		}

		if(cp > 1e-7)
		{
			if(matrix_confidence_.col(index_way_last_mapped).maxCoeff() > 0.1)
			{
				num_q1 += 1;
				q1 += cp;
			}
			else
				num_add += 1;
		}
		else
		{
			if(matrix_confidence_.col(index_way_last_mapped).maxCoeff() < 0.1)
				num_add += 1;
		}
	}
	if(mapped.size() - num_add > 1e-7)
		q1 = q1 / (mapped.size()- num_add);
	else
		q1 = matrix_confidence_(index_way_main, index_way_mapped);

	//q2
	// j,i,j,a
	std::vector<PWay>::iterator it_last_main = main.begin();
	for(; it_last_main != main.end(); ++it_last_main)
	{
		PWay way_last_main = *it_last_main;
		std::map<std::string, int>::iterator it_find_last_main = ids_main.find(way_last_main->osm_id);
		if (it_find_last_main == ids_main.end())
			continue;
		int index_way_last_main = ids_main[way_last_main->osm_id];
		double confidence_j = matrix_confidence_(index_way_last_main, index_way_mapped);
		if(confidence_j > 1e-7 )
		{
			double d1(0.);
			double d2(0.);
			double rate(0.);
			std::string key(way_last_main->osm_id + "$" + way_mapped->osm_id);
			MapWayMatchingPair::iterator it_find = map_way_matching_pair_.find(key);
			if (it_find != map_way_matching_pair_.end())
			{
				d1 = it_find->second->dis_mapping_sta;
				d2 = it_find->second->dis_mapping_end;
				if(it_find->second->length_mapping_ln1 < 1e-7)
					continue;
				rate = fabs(1.0 - fabs(it_find->second->length_mapping_ln1 - it_find->second->length_mapping_ln2) / std::max(it_find->second->length_mapping_ln1, it_find->second->length_mapping_ln2));
			}
			else
				continue;
			double q_dis = 0.0;
			double q_dir = 0.0;
			if((dis_hausdorff0 + dis_hausdorff1) > 1e-7)
			{
				q_dis = fabs((d1 + d2) / 2 - (dis_hausdorff0 + dis_hausdorff1) / 2);
			}
			double d_dir(0.);
//			if(last)
//				d_dir = fabs(GeosTool::calculate_include_angle(0, GeosTool::calculate_include_angle(way_last_main->end_angle, angle_start_main)));
//			else
//				d_dir = fabs(GeosTool::calculate_include_angle(0, GeosTool::calculate_include_angle(way_last_main->sta_angle, angle_end_main)));
			q_dir = d_dir;

			double a = -0.2*std::log(2);
			double b = 4 * std::log(2);
			double dis = 1.0 / (1.0 + std::pow(M_E, -(a * q_dis + b)));
			a = -8*std::log(2);
			b = 6 * std::log(2);
			double dir = 1.0 / (1.0 + std::pow(M_E, -(a * q_dir + b)));

			if (q_dis < 1e-7)
				dis = 1.0;
			if (q_dir < 1e-7)
				dir = 1.0;

			double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * rate;

			if(debug)
			{
				std::cout << "---------------------i,j,i,b-------------------------" << std::endl;
				std::cout << "		main_last_way: " << way_last_main->osm_id << ",	mapped_last_way: " << way_mapped->osm_id << std::endl;
				std::cout << "		main_last_way_sta_angle: " << way_last_main->sta_angle << ",	main_last_way_end_angle: " << way_last_main->end_angle << std::endl;
				std::cout << "		mapped_last_way_sta_angle: " << way_mapped->sta_angle << ",	mapped_last_way_end_angle: " << way_mapped->end_angle << std::endl;
				std::cout << "		dis1: " << d1 << ",	dis2: " << d2 << std::endl;
				std::cout << "		d_dir: " << d_dir << std::endl;
				std::cout << "		q_dis: " << q_dis << ",	q_dir: " << q_dir << std::endl;
				std::cout << "		length_mapping_ln2: " << it_find->second->length_mapping_ln2 << ",	length_mapping_ln1: " << it_find->second->length_mapping_ln1 << std::endl;
				std::cout << "		rate:  " << rate  << std::endl;
				std::cout << "		c: " << c << ",	cp(q2): " << c * confidence_j << std::endl;
			}

			q2 = std::max(q2, c * confidence_j);
		}
	}
//	if(q1 > 1e-7 || q2 > 1e-7)
//	{
//		if(q1 > 1e-7 && q2 > 1e-7)
//		{
//			q = (q1 + q2) / (main.size() - num_add + 1);
//		}
//		else if(q1 > 1e-7)
//		{
//			if(main.size() - num_add > 1e-7)
//			{
//				q = q1 / (main.size() - num_add);
//			}
//			else
//			{
//				q = matrix_confidence_(index_way_main, index_way_mapped);
//			}
//		}
//		else
//		{
//			q = q2;
//		}
//	}
//	else
//		q = 0.0;
	q = std::max(q1, q2);
	if(debug)
	{
		std::cout << "JI: num_q1: " << num_q1 << std::endl;
		std::cout << "q1: " << q1 << ", q2: " << std::endl;
		std::cout << "\n	JI Result:  " << q << ",	num_added: " << num_add << std::endl;
	}
}

void MatchingProbabilityCalculator::initEigen_support(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
{
	bool debug(false);
	MapWayMatchingPair::iterator it_wmp = map_way_matching_pair_.begin();
	for (; it_wmp != map_way_matching_pair_.end(); ++it_wmp)
	{
		PWayMatchingPair pWayMatchingPair = it_wmp->second;
		PWay way_main = pWayMatchingPair->main_way;
		PWay way_mapped = pWayMatchingPair->mapped_way;

		int index_way_main = ids_main[way_main->osm_id];
		int index_way_mapped = ids_mapped[way_mapped->osm_id];

		double q_last_ij(0.);
		double q_next_ij(0.);
		double q_last_ji(0.);
		double q_next_ji(0.);
		MatchingProbabilityCalculator::calQIJ(q_last_ij, pWayMatchingPair, ids_main, ids_mapped, true);
		MatchingProbabilityCalculator::calQIJ(q_next_ij, pWayMatchingPair, ids_main, ids_mapped, false);
		MatchingProbabilityCalculator::calQJI(q_last_ji, pWayMatchingPair, ids_main, ids_mapped, true);
		MatchingProbabilityCalculator::calQJI(q_next_ji, pWayMatchingPair, ids_main, ids_mapped, false);
		if(debug)
		{
			std::cout << "		==================calQ================" << std::endl;
			std::cout << "			q_last_ij: " << q_last_ij <<  "  q_next_ij: " << q_next_ij << "  q_last_ji: " << q_last_ji << "  q_next_ji: " << q_next_ji << std::endl;
		}
		double q_ij = q_last_ij + q_next_ij;
		double q_ji = q_last_ji + q_next_ji;

		double q = sqrt(q_ij * q_ji);

		matrix_support_(index_way_main, index_way_mapped) = q;
	}
}

int MatchingProbabilityCalculator::updateEigen_confidence(std::string fname, std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
{
	int num_change(0);
	std::ofstream f_out(fname);
	MapWayMatchingPair::iterator it_wmp = map_way_matching_pair_.begin();
	for (; it_wmp != map_way_matching_pair_.end(); ++it_wmp)
	{
		PWayMatchingPair pWayMatchingPair = it_wmp->second;

		PWay way_main = pWayMatchingPair->main_way;
		PWay way_mapped = pWayMatchingPair->mapped_way;

		int index_way_main = ids_main[way_main->osm_id];
		int index_way_mapped = ids_mapped[way_mapped->osm_id];

		double p_old = matrix_confidence_(index_way_main, index_way_mapped);

 		if(pWayMatchingPair->confidence_type != MATCH_XGBOOST && pWayMatchingPair->confidence_type != DISBELIEF)
		{
			f_out << way_main->osm_id + "\t" + way_mapped->osm_id + "\t" + std::to_string(p_old) + "\n";
			continue;
		}

		// if(p_old < 1e-7 || way_main->vec_last_way.size() < 1 || way_main->vec_next_way.size() < 1 || way_mapped->vec_last_way.size() < 1 || way_mapped->vec_next_way.size() < 1)
		if(p_old < 1e-7)
		{
			f_out << way_main->osm_id + "\t" + way_mapped->osm_id + "\t" + std::to_string(p_old) + "\n";
			continue;
		}
		double q_old = matrix_support_(index_way_main, index_way_mapped);
		double p_i_new = (p_old + q_old) / (1 + matrix_support_.row(index_way_main).sum());
		double p_j_new = (p_old + q_old) / (1 + matrix_support_.col(index_way_mapped).sum());
		double p_new = 2 * p_i_new * p_j_new / (p_i_new + p_j_new);
		matrix_confidence_(index_way_main, index_way_mapped) = p_new;
		if(fabs(p_new-p_old) > 0.1)
		{
			num_change += 1;
//			std::cout << "		main_id: " << way_main->osm_id << ",	mapped_id: " << way_mapped->osm_id << std::endl;
//			std::cout << "			p_new: " << p_new << ",	p_old: " << p_old << ",	change: " << p_new - p_old << std::endl;
		}
//		if (way_main->osm_id == "w1712511798765" && way_mapped->osm_id == "w12034303")
//		{
//			std::cout << "			p_new: " << p_new << ",	p_old: " << p_old << ",	change: " << p_new - p_old << std::endl;
//		}

		f_out << way_main->osm_id + "\t" + way_mapped->osm_id + "\t" + std::to_string(p_old) + "\t" + std::to_string(p_new) << "\t" << std::to_string(p_new - p_old) + "\n";
	}
	f_out.close();
	return num_change;
}

int MatchingProbabilityCalculator::updateEigen_support(std::string fname, std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
{
	int num_change(0);
	std::ofstream f_out(fname);
	MapWayMatchingPair::iterator it_wmp = map_way_matching_pair_.begin();
	for (; it_wmp != map_way_matching_pair_.end(); ++it_wmp)
	{
		PWayMatchingPair pWayMatchingPair = it_wmp->second;
		PWay way_main = pWayMatchingPair->main_way;
		PWay way_mapped = pWayMatchingPair->mapped_way;

		int index_way_main = ids_main[way_main->osm_id];
		int index_way_mapped = ids_mapped[way_mapped->osm_id];
		double q_old = matrix_support_(index_way_main, index_way_mapped);
		if(q_old < 1e-7)
		{
			continue;
		}

		double q_last_ij(0.);
		double q_next_ij(0.);
		double q_last_ji(0.);
		double q_next_ji(0.);
		MatchingProbabilityCalculator::calQIJ(q_last_ij, pWayMatchingPair, ids_main, ids_mapped, true);
		MatchingProbabilityCalculator::calQIJ(q_next_ij, pWayMatchingPair, ids_main, ids_mapped, false);
		MatchingProbabilityCalculator::calQJI(q_last_ji, pWayMatchingPair, ids_main, ids_mapped, true);
		MatchingProbabilityCalculator::calQJI(q_next_ji, pWayMatchingPair, ids_main, ids_mapped, false);
		double q_ij = q_last_ij + q_next_ij;
		double q_ji = q_last_ji + q_next_ji;

		double q_new = sqrt(q_ij * q_ji);
		matrix_support_(index_way_main, index_way_mapped) = q_new;
//		if(q_new != q_old)
		if(fabs(q_new-q_old) > 0.1)
		{
			num_change += 1;
		}
		f_out << way_main->osm_id + "\t" + way_mapped->osm_id + "\t" + std::to_string(q_old) + "\t" + std::to_string(q_new) << "\t" << std::to_string(q_new - q_old) + "\n";
	}

	f_out.close();
	return num_change;
}

void MatchingProbabilityCalculator::update()
{
	//	/Users/didi/Documents/work/data/result_0519/iteration
	bool debug(true);
	if(debug)
	{
		std::cout << "------------------initEigen--------------------" << std::endl;
		std::cout << "	ids:	ids_main_: " << ids_main_.size() << ", ids_mapped_: " << ids_mapped_.size() << std::endl;
		std::cout << "	SIZE matrix_confidence_:	row:" << matrix_confidence_.rows() << ", cols" << matrix_confidence_.cols() << ", sum: " << matrix_confidence_.sum() << std::endl;
//		MatchingProbabilityCalculator::outMatrix("/Users/didi/Documents/work/data/result_0519/iteration/p0.txt", matrix_confidence_);
//		MatchingProbabilityCalculator::outIds("/Users/didi/Documents/work/data/result_0519/iteration/ids_main_.txt", ids_main_);
//		MatchingProbabilityCalculator::outIds("/Users/didi/Documents/work/data/result_0519/iteration/ids_mapped_.txt", ids_mapped_);
	}
	MatchingProbabilityCalculator::initEigen_support(ids_main_, ids_mapped_);
	if(debug)
	{
		std::cout << "------------------initEigen_support--------------------" << std::endl;
		std::cout << "	SIZE matrix_support_:	row:" << matrix_support_.rows() << ", cols" << matrix_support_.cols() << ", sum: " << matrix_support_.sum()<< std::endl;
//		MatchingProbabilityCalculator::outMatrix("/Users/didi/Documents/work/data/result_0519/iteration/q0.txt", matrix_support_);
//		MatchingProbabilityCalculator::outPairConfidence("/Users/didi/Documents/work/data/result_0519/iteration/Pair0.txt", ids_main_, ids_mapped_);
	}
	if(debug)
		std::cout << "------------------start update--------------------" << std::endl;
	int i = 0;
	int num_change(INT_MAX);
//	while(num_change > CHANGED)
	while(num_change > 0)
	{
		if(debug)
		{
			std::cout << "	N: " << ++i << "	num_change:  " << num_change << std::endl;
//			int index_way_main = ids_main_["w1712508846153"];
//			int index_way_mapped = ids_mapped_["w1709007230079"];
//			std::cout << "		=======confidence: " << matrix_confidence_(index_way_main, index_way_mapped) << std::endl;
		}

		std::string fname1 = "/Users/didi/Documents/work/data/result_0519/iteration/p"+ std::to_string(i) +".txt";
		int n_change_confidence = MatchingProbabilityCalculator::updateEigen_confidence(fname1, ids_main_, ids_mapped_);
		std::cout << "	------------------updateEigen_confidence--------------------" << std::endl;

		std::string fname2 = "/Users/didi/Documents/work/data/result_0519/iteration/q"+ std::to_string(i) +".txt";
		int n_change_support = MatchingProbabilityCalculator::updateEigen_support(fname2, ids_main_, ids_mapped_);
		std::cout << "	------------------updateEigen_support--------------------" << std::endl;

//		num_change = std::max(n_change_confidence, n_change_support);
		if (n_change_confidence >= num_change)
			break;

		num_change = n_change_confidence;
		if(debug)
		{
//			MatchingProbabilityCalculator::outMatrix("/Users/didi/Documents/work/data/result_0519/iteration/p"+ std::to_string(i) +".txt", matrix_confidence_);
//			MatchingProbabilityCalculator::outMatrix("/Users/didi/Documents/work/data/result_0519/iteration/q"+ std::to_string(i) +".txt", matrix_support_);
//			MatchingProbabilityCalculator::outPairConfidence("/Users/didi/Documents/work/data/result_0519/iteration/Pair"+ std::to_string(i) +".txt", ids_main, ids_mapped);
		}
		return;
	}
	std::cout << "OVER!" << std::endl;
}

void MatchingProbabilityCalculator::outMatrix(std::string fname, Eigen::MatrixXd datas)
{
	std::ofstream f_p(fname);
	for(int i=0; i<datas.rows(); i++)
	{
		Eigen::MatrixXd row = datas.row(i);
		std::string out("");
		for(int j=0; j<row.size(); j++)
		{
			out += std::to_string(row(0, j)) + "\t";
		}
		out += "\n";
		f_p << out;
	}
	f_p.close();
}

void MatchingProbabilityCalculator::outIds(std::string fname, std::map<std::string, int>& ids)
{
	std::ofstream f_id(fname);
	std::map<std::string, int>::iterator id = ids.begin();
	for (; id != ids.end(); ++id)
	{
		f_id << id->first + "\t" + std::to_string(id->second) + "\n";
	}
	f_id.close();
}

void MatchingProbabilityCalculator::outPairConfidence(std::string fname, std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
{
	std::ofstream f_p(fname);
	MapWayMatchingPair::iterator it_wmp = map_way_matching_pair_.begin();
	for (; it_wmp != map_way_matching_pair_.end(); ++it_wmp)
	{
		PWayMatchingPair pWayMatchingPair = it_wmp->second;
		PWay way_main = pWayMatchingPair->main_way;
		PWay way_mapped = pWayMatchingPair->mapped_way;
		int index_main = ids_main[way_main->osm_id];
		int index_mapped = ids_mapped[way_mapped->osm_id];
		f_p << std::to_string(matrix_confidence_(index_main, index_mapped)) + "\t" + std::to_string(matrix_support_(index_main, index_mapped)) + "\n";
	}
	f_p.close();
}


