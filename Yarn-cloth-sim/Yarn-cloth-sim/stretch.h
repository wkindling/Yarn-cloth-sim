#pragma once
#ifndef STRETCH_H
#define STRETCH_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "node.h"

class StretchSpring
{
public:
	StretchSpring(Node* n0, Node* n1, double _Y, double _R, YarnType type);
	virtual ~StretchSpring();

	void solve();

	void solveU();
	void solveV();

	void draw();

public:
	Node *node0, *node1;

	double stretchEnergy;
	double Y, R, Ks;

	YarnType springType;
};


#endif