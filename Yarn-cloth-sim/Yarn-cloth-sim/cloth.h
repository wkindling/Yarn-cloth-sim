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
	Cloth(int w, int h, double _R, double _L, double _mu, double _rho, double _Y, double _B, double _S, double _Kc, double _Kf);
	virtual ~Cloth();

	void build();
	void step(double h);
	void draw();

	void computeForce(Eigen::Vector3d gravity);
	void computeInertia();
	void solve(double h);
	
public:
	int width, height;
	double R, L, mu, rho;
	double Y, B, S, Kc, Kf;

	int DoF;

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