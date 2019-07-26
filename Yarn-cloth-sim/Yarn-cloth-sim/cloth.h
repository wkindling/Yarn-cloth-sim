#pragma once
#ifndef CLOTH_H
#define CLOTH_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "node.h"

class ShearSpring;
class StretchSpring;
class ParallelContactSpring;
class BendSpring;

class Cloth
{
public:
	Cloth(int w, int h, double _R, double _L);
	virtual ~Cloth();

	void build();
	void step();
	void draw();

	void computeForce();
	void computeInertia();
	void solve();
	
public:
	int width, height;
	double R, L, mu, rho;
	double Y, B, S, Kc, Kf;

	std::vector<Node*> nodes;

	std::vector<ShearSpring*> shear_springs;
	std::vector<StretchSpring*> stretch_springs;
	std::vector<ParallelContactSpring*> parallel_contact_springs;
	std::vector<BendSpring*> bend_springs;


	Eigen::VectorXd v;
	Eigen::VectorXd f;
	Eigen::SparseMatrix<double> M;
	Eigen::SparseMatrix<double> K;


};



#endif