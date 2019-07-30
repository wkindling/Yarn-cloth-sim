#pragma once
#ifndef STRETCH_H
#define STRETCH_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#include "node.h"

class StretchSpring
{
public:
	StretchSpring(Node* n0, Node* n1, double _Y, double _R, YarnType type);
	virtual ~StretchSpring();

	void solve(std::vector<T>& _K, Eigen::VectorXd& f, int nodes_size);

	void solveU(std::vector<T>& _K, Eigen::VectorXd& f, int nodes_size);
	void solveV(std::vector<T>& _K, Eigen::VectorXd& f, int nodes_size);

	void draw();

public:
	Node *node0, *node1;

	double stretchEnergy;
	double Y, R, Ks;

	YarnType springType;
};


#endif