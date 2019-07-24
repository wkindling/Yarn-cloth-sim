#pragma once
#ifndef BEND_H
#define BEND_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "node.h"

class BendSpring
{
public:
	BendSpring(Node* n0, Node* n1, Node* n2, double _B, double _R, YarnType type);
	virtual ~BendSpring();

	void solve(std::vector<T>& _K, Eigen::VectorXd& f);

	void solveU(std::vector<T>& _K, Eigen::VectorXd& f);
	void solveV(std::vector<T>& _K, Eigen::VectorXd& f);

public:
	Node *node0, *node1, *node2;
	
	double bendEnergy;
	double B, R, Kb;

	YarnType springType;
};

#endif