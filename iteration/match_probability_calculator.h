/*
 * matching_probability_calculator.h
 *
 *  Created on: 2018年5月10日
 *      Author: liping
 */

#ifndef SRC_ITERATION_MATCHING_PROBABILITY_CALCULATOR_H_
#define SRC_ITERATION_MATCHING_PROBABILITY_CALCULATOR_H_

#include "feature_probability_calculator.h"

class MatchingProbabilityCalculator
{
public:
	MatchingProbabilityCalculator(WayOperator* way_operator1, WayOperator* way_operator2, const MapWayMatchingPair& map_way_matching_pair);
	virtual ~MatchingProbabilityCalculator();

	void calculator();

private:
	void init();
	void initInitialProbability(); 			// 初始化初始概率
	void calculatorAdjacencyProbability(); 	// 计算支持系数
	void calculatorMatchingProbability(); 	// 计算匹配概率
	void updateMatchingPairProbability();	//

private:
	WayOperator* way_operator1_;
	WayOperator* way_operator2_;
	MapWayMatchingPair map_way_matching_pair_;
	Eigen::MatrixXd* matrix_initial_probability_;	// 初始概率
	Eigen::MatrixXd* matrix_adjacency_probability_; 	// 兼容系数和支持系数
	Eigen::MatrixXd* matrix_matched_probability_; 	// 匹配概率
	FeatureProbabilityCalculator* feature_probability_calculator_;
};

#endif /* SRC_ITERATION_MATCHING_PROBABILITY_CALCULATOR_H_ */
