/*
 * feature_probability_calculator.h
 *
 *  Created on: 2018年5月10日
 *      Author: yinweijun
 */

#ifndef SRC_ITERATION_FEATURE_PROBABILITY_CALCULATOR_H_
#define SRC_ITERATION_FEATURE_PROBABILITY_CALCULATOR_H_

#include "../road/data_operator.h"
#include "../matching/way_matching_operator.h"
#include <vector>
#include <math.h>
#include "Eigen/Dense"

//struct Hausdorff
//{
////	double length;
//	double distance;
//	double direction;
//	Hausdorff() : /*length(-1.), */distance(-1.), direction(-1.) {}
//};
//
//struct Feature
//{
//	double length;
//	double distance;
//	double direction;
//	Hausdorff hausdorff[2];// 0:i->j 1:j->i
//	Feature() : length(-1.), distance(-1.), direction(-1.) {}
//
//};

class FeatureProbabilityCalculator
{
public:
	FeatureProbabilityCalculator(WayOperator* way_operator1, WayOperator* way_operator2, const MapWayMatchingPair& map_way_matching_pair);
	virtual ~FeatureProbabilityCalculator();

	void calculator();

private:
	void init();
	void initEigen(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped);			// 初始化概率矩阵
	void calQ(double& q, PWayMatchingPair& pWayMatchingPair, std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped, bool i_j, bool last);		// 计算匹配对的支撑系数，last = true，前继， last=false，后继
	void initEigen_support(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped);
	void updateEigen_confidence(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped);
	void updateEigen_support(std::map<std::string, int>& ids_main, std::map<std::string, int>& ids_mapped);
	void update(); 			// 计算特征矩阵

private:
	Eigen::MatrixXd matrix_confidence_;
	Eigen::MatrixXd matrix_support_;
	Eigen::MatrixXd matrix_hausdorff_dis0_;
	Eigen::MatrixXd matrix_hausdorff_dis1_;
	Eigen::MatrixXd matrix_hausdorff_dir0_;
	Eigen::MatrixXd matrix_hausdorff_dir1_;

private:
	WayOperator* way_operator1_;
	WayOperator* way_operator2_;
	MapWayMatchingPair map_way_matching_pair_;
};

#endif /* SRC_ITERATION_FEATURE_PROBABILITY_CALCULATOR_H_ */
