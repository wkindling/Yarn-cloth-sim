#pragma once
#ifndef PARALLEL_CONTACT_H
#define PARALLEL_CONTACT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#include "node.h"

class ParallelContactSpring
{
public:
	ParallelContactSpring(Node* n0, Node* n1, double _Kc, double _R, double _L, YarnType type);
	virtual ~ParallelContactSpring();

	void solve(std::vector<T>& _K, Eigen::VectorXd& f);
	void solveU(std::vector<T>& _K, Eigen::VectorXd& f);
	void solveV(std::vector<T>& _K, Eigen::VectorXd& f);

public:
	Node *node0, *node1;
	double Kc;
	double L, R;
	double d;

	double parallelContactEnergy;

	YarnType springType;
};

#endif