#pragma once
#ifndef SHEAR_H
#define SHEAR_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "node.h"

class ShearSpring
{
public:
	ShearSpring(Node* n0, Node* n1, Node* n3, double _S, double _R, double _L);
	virtual ~ShearSpring();

	void solve();

public:
	Node *node0, *node1, *node3;

	double shearEnergy;
	double S, R, L, Kx;



};



#endif