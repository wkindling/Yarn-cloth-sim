#pragma once
#ifndef SHEAR_H
#define SHEAR_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#include "node.h"

class ShearSpring
{
public:
	ShearSpring(Node* n0, Node* n1, Node* n3, double _S, double _R, double _L);
	virtual ~ShearSpring();

	void solve(std::vector<T>&_K, Eigen::VectorXd& f, double h);

public:
	Node *node0, *node1, *node3;

	double shearEnergy;
	double S, R, L, Kx;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};



#endif