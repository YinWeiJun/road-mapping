/*
 * feature_probability_calculator.cpp
 *
 *  Created on: 2018年5月10日
 *      Author: yinweijun
 */

#include "feature_probability_calculator.h"

FeatureProbabilityCalculator::FeatureProbabilityCalculator(
		WayOperator* way_operator1, WayOperator* way_operator2,
		const MapWayMatchingPair& map_way_matching_pair)
:  way_operator1_(way_operator1)
, way_operator2_(way_operator2)
, map_way_matching_pair_(map_way_matching_pair)
{
	init();
}

FeatureProbabilityCalculator::~FeatureProbabilityCalculator()
{
}

void FeatureProbabilityCalculator::calculator()
{
	update();
}

void FeatureProbabilityCalculator::init()
{
	int rows = way_operator1_->map_way_osm_id_.size();
	int cols = way_operator2_->map_way_osm_id_.size();
	matrix_confidence_.Zero(rows, cols);
	matrix_support_.Zero(rows, cols);
	matrix_hausdorff_dis0_.Zero(rows, cols);
	matrix_hausdorff_dis1_.Zero(rows, cols);
	matrix_hausdorff_dir0_.Zero(rows, cols);
	matrix_hausdorff_dir1_.Zero(rows, cols);
}

void FeatureProbabilityCalculator::initEigen(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
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
			index_mapped = ids_main[way_mapped->osm_id];
		}
		else
		{
			ids_mapped.insert(make_pair(way_mapped->osm_id, num_cols));
			num_cols += 1;
		}

		matrix_confidence_(index_main, index_mapped) = pWayMatchingPair->confidence;
		matrix_hausdorff_dis0_(index_main, index_mapped) = pWayMatchingPair->dis_hausdorff_ln1;
		matrix_hausdorff_dis1_(index_main, index_mapped) = pWayMatchingPair->dis_hausdorff_ln2;
		matrix_hausdorff_dir0_(index_main, index_mapped) = pWayMatchingPair->angle_hausdorff_ln1;
		matrix_hausdorff_dir1_(index_main, index_mapped) = pWayMatchingPair->angle_hausdorff_ln2;
	}
}

void FeatureProbabilityCalculator::calQ(double& q, PWayMatchingPair& pWayMatchingPair, std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped, bool i_j, bool last)
{
	PWay way_main = pWayMatchingPair->main_way;
	PWay way_mapped = pWayMatchingPair->mapped_way;

	if(!i_j)
	{
		way_main = pWayMatchingPair->mapped_way;
		way_mapped = pWayMatchingPair->main_way;
	}
	int index_way_main = ids_main[way_main->osm_id];
	int index_way_mapped = ids_mapped[way_mapped->osm_id];

	double dis_hausdorff0 = matrix_hausdorff_dis0_(index_way_main, index_way_mapped);
	double dis_hausdorff1 = matrix_hausdorff_dis1_(index_way_main, index_way_mapped);

	double dir_hausdorff0 = matrix_hausdorff_dir0_(index_way_main, index_way_mapped);
	double dir_hausdorff1 = matrix_hausdorff_dir1_(index_way_main, index_way_mapped);

	double angle_start_main = way_main->sta_angle;
	double angle_end_main = way_main->end_angle;
	double angle_start_mapped = way_mapped->sta_angle;
	double angle_end_mapped = way_mapped->end_angle;

	std::vector<PWay> main;
	std::vector<PWay> mapped;

	if(last)
	{
		main = way_main->vec_last_way;
		mapped = way_mapped->vec_last_way;
	}
	else
	{
		main = way_main->vec_next_way;
		mapped = way_mapped->vec_next_way;
	}

	double q1(0.);
	double q2(0.);

	// q1
	int num_q1(0);
	std::vector<PWay>::iterator it_main = main.begin();
	for(; it_main != main.end(); ++it_main)
	{
		PWay way_last_main = *it_main;
		int index_way_last_main = ids_main[way_last_main->osm_id];
		std::vector<PWay>::iterator it_last_mapped = mapped.begin();
		double cp(0.);
		for(; it_last_mapped != mapped.end(); ++it_last_mapped)
		{
			PWay way_last_mapped = *it_last_mapped;
			int index_way_last_mapped = ids_mapped[way_last_mapped->osm_id];
			// i,j,a,b
			double confidence = matrix_confidence_(index_way_last_main, index_way_last_mapped);
			if(confidence > 1e-7 )
			{
				double dis_hausdorff_last0 = matrix_hausdorff_dis0_(index_way_last_main, index_way_last_mapped);
				double dis_hausdorff_last1 = matrix_hausdorff_dis1_(index_way_last_main, index_way_last_mapped);

				double q_dis = (dis_hausdorff_last0 + dis_hausdorff_last1) / (dis_hausdorff0 + dis_hausdorff1);

				double d_dir(0.);
				if(last)
					d_dir = fabs((way_last_main->end_angle - angle_start_main) - (way_last_mapped->end_angle - angle_start_mapped));
				else
					d_dir = fabs((way_last_main->sta_angle - angle_end_main) - (way_last_mapped->sta_angle - angle_end_mapped));
				double q_dir = 2 * d_dir / (dir_hausdorff0 + dir_hausdorff1);

				double dis = 1.0 / (1.0 + q_dis * q_dis);
				double dir = 1.0 / (1.0 + q_dir * q_dir);

				double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * (way_last_main->length / way_last_mapped->length);

				cp = std::max(cp, c * confidence);
			}
		}
		// i,j,a.j
		double confidence_j = matrix_confidence_(index_way_last_main, index_way_mapped);
		if(confidence_j > 1e-7 )
		{
			double dis_hausdorff_last0 = matrix_hausdorff_dis0_(index_way_last_main, index_way_mapped);
			double dis_hausdorff_last1 = matrix_hausdorff_dis1_(index_way_last_main, index_way_mapped);

			double q_dis = (dis_hausdorff_last0 + dis_hausdorff_last1) / (dis_hausdorff0 + dis_hausdorff1);

			double d_dir(0.);
			if(last)
				d_dir = fabs((way_last_main->end_angle - angle_start_main) - 0);
			else
				d_dir = fabs((way_last_main->sta_angle - angle_end_main) - 0);
			double q_dir = 2 * d_dir / (dir_hausdorff0 + dir_hausdorff1);

			double dis = 1.0 / (1.0 + q_dis * q_dis);
			double dir = 1.0 / (1.0 + q_dir * q_dir);

			double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * (way_last_main->length / way_mapped->length);

			cp = std::max(cp, c * confidence_j);
		}

		if(cp > 1e-7)
		{
			num_q1 += 1;
			q1 += cp;
		}
	}
	q1 = q1 / num_q1;

	//q2
	// i,j,i,b
	std::vector<PWay>::iterator it_last_mapped = mapped.begin();
	for(; it_last_mapped != mapped.end(); ++it_last_mapped)
	{
		PWay way_last_mapped = *it_last_mapped;
		int index_way_last_mapped = ids_mapped[way_last_mapped->osm_id];
		double confidence_i = matrix_confidence_(index_way_main, index_way_last_mapped);
		if(confidence_i > 1e-7 )
		{
			double dis_hausdorff_last0 = matrix_hausdorff_dis0_(index_way_main, index_way_last_mapped);
			double dis_hausdorff_last1 = matrix_hausdorff_dis1_(index_way_main, index_way_last_mapped);

			double q_dis = (dis_hausdorff_last0 + dis_hausdorff_last1) / (dis_hausdorff0 + dis_hausdorff1);

			double d_dir(0.);
			if(last)
				d_dir = fabs(0 - (way_last_mapped->end_angle - angle_start_mapped));
			else
				d_dir = fabs(0 - (way_last_mapped->sta_angle - angle_end_mapped));
			double q_dir = 2 * d_dir / (dir_hausdorff0 + dir_hausdorff1);

			double dis = 1.0 / (1.0 + q_dis * q_dis);
			double dir = 1.0 / (1.0 + q_dir * q_dir);

			double c = (2.0 / ((1.0 / dis) + (1.0 / dir))) * (way_main->length / way_last_mapped->length);

			q2 = std::max(q2, c * confidence_i);
		}
	}

	q = std::max(q1, q2);
}

void FeatureProbabilityCalculator::initEigen_support(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
{
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
		FeatureProbabilityCalculator::calQ(q_last_ij, pWayMatchingPair, ids_main, ids_mapped, true, true);
		FeatureProbabilityCalculator::calQ(q_next_ij, pWayMatchingPair, ids_main, ids_mapped, true, false);
		FeatureProbabilityCalculator::calQ(q_last_ji, pWayMatchingPair, ids_main, ids_mapped, false, true);
		FeatureProbabilityCalculator::calQ(q_next_ji, pWayMatchingPair, ids_main, ids_mapped, false, false);
		double q_ij = q_last_ij + q_next_ij;
		double q_ji = q_last_ji + q_next_ji;

		double q = sqrt(q_ij * q_ji);
		matrix_confidence_(index_way_main, index_way_mapped) = q;
	}
}

void FeatureProbabilityCalculator::updateEigen_confidence(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
{
	MapWayMatchingPair::iterator it_wmp = map_way_matching_pair_.begin();
	for (; it_wmp != map_way_matching_pair_.end(); ++it_wmp)
	{
		PWayMatchingPair pWayMatchingPair = it_wmp->second;
		PWay way_main = pWayMatchingPair->main_way;
		PWay way_mapped = pWayMatchingPair->mapped_way;

		int index_way_main = ids_main[way_main->osm_id];
		int index_way_mapped = ids_mapped[way_mapped->osm_id];

		double p_old = matrix_confidence_(index_way_main, index_way_mapped);
		double q_old = matrix_support_(index_way_main, index_way_mapped);

		double p_i_new = (p_old + q_old) / (1 + matrix_support_.row(index_way_main).sum());
		double p_j_new = (p_old + q_old) / (1 + matrix_support_.col(index_way_mapped).sum());

		double p_new = 2 * p_i_new * p_j_new / (p_i_new + p_j_new);

		matrix_confidence_(index_way_main, index_way_mapped) = p_new;
	}
}

void FeatureProbabilityCalculator::updateEigen_support(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped)
{
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
		FeatureProbabilityCalculator::calQ(q_last_ij, pWayMatchingPair, ids_main, ids_mapped, true, true);
		FeatureProbabilityCalculator::calQ(q_next_ij, pWayMatchingPair, ids_main, ids_mapped, true, false);
		FeatureProbabilityCalculator::calQ(q_last_ji, pWayMatchingPair, ids_main, ids_mapped, false, true);
		FeatureProbabilityCalculator::calQ(q_next_ji, pWayMatchingPair, ids_main, ids_mapped, false, false);
		double q_ij = q_last_ij + q_next_ij;
		double q_ji = q_last_ji + q_next_ji;

		double q_new = sqrt(q_ij * q_ji);
		matrix_support_(index_way_main, index_way_mapped) = q_new;
	}
}

void FeatureProbabilityCalculator::update()
{
	std::map<std::string, int> ids_main;
	std::map<std::string, int> ids_mapped;

	FeatureProbabilityCalculator::initEigen(ids_main, ids_mapped);
	FeatureProbabilityCalculator::initEigen_support(ids_main, ids_mapped);

	Eigen::MatrixXd tmp;
	tmp.Zero(ids_main.size(), ids_mapped.size());
	Eigen::MatrixXd old = matrix_confidence_;
	Eigen::MatrixXd diff = matrix_confidence_;
	while(diff != tmp)
	{
		FeatureProbabilityCalculator::updateEigen_confidence(ids_main, ids_mapped);
		FeatureProbabilityCalculator::updateEigen_support(ids_main, ids_mapped);
		diff = matrix_confidence_ - old;
		old = matrix_confidence_;
	}
}
