#pragma once
#ifndef CLOTH_H
#define CLOTH_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>

class Node;
class ShearSpring;
class StretchSpring;
class ParallelContactSpring;
class BendSpring;

class Cloth
{
public:
	Cloth(int w, int h, double _R, double _L);
	virtual ~Cloth();

	void draw();

public:
	int width, height;
	double R, L;
	double Y, B, S, Kc;

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