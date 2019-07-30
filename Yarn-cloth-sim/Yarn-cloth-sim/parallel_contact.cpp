#include "parallel_contact.h"
#include <iostream>

using namespace std;
using namespace Eigen;

ParallelContactSpring::ParallelContactSpring(Node* n0, Node* n1, double _Kc, double _R, double _L, YarnType type)
{
	//To ensure that 0->1 is the positive direction of the rod
	if (type == Warp)
	{
		node0 = n0->u < n1->u ? n0 : n1;
		node1 = n1->u > n0->u ? n1 : n0;
	}
	else if (type == Weft)
	{
		node0 = n0->v < n1->v ? n0 : n1;
		node1 = n1->v > n0->v ? n1 : n0;
	}

	Kc = _Kc;
	L = _L;
	R = _R;

	d = 4 * R;
	parallelContactEnergy = 0;

	springType = type;
}

ParallelContactSpring::~ParallelContactSpring() {}

void ParallelContactSpring::solve(vector<T>& _K, VectorXd& f, int nodes_size)
{
	if (springType == Weft) solveV(_K, f, nodes_size);

	else if (springType == Warp) solveU(_K, f, nodes_size);
	

	return;
}

void ParallelContactSpring::solveU(vector<T>& _K, VectorXd& f, int nodes_size)
{
	double l = (node1->position - node0->position).norm();
	double delta_u = abs(node1->u - node0->u);

	if (delta_u >= d) return;

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d w = (node1->position - node0->position) / delta_u;

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d1 * d1.transpose();

	parallelContactEnergy = 0.5*Kc*L*(delta_u - d)*(delta_u - d);

	if (!node0->onBorder)
	{
		int cross_index0 = nodes_size * 3 + node0->cross_index * 2;
		double Fu0 = Kc * L*(delta_u - d);
		f(cross_index0) += Fu0;

		double Fu0du0 = -Kc * L;
		_K.push_back(T(cross_index0, cross_index0, Fu0du0));
	}

	if (!node1->onBorder)
	{
		int cross_index1 = nodes_size * 3 + node1->cross_index * 2;
		double Fu1 = -Kc * L*(delta_u - d);
		f(cross_index1) += Fu1;

		double Fu1du1 = -Kc * L;
		_K.push_back(T(cross_index1, cross_index1, Fu1du1));
	}

	if (!node0->onBorder && !node1->onBorder)
	{
		int cross_index0 = nodes_size * 3 + node0->cross_index * 2;
		int cross_index1 = nodes_size * 3 + node1->cross_index * 2;

		double Fu0du1 = Kc * L;
		double Fu1du0 = Fu0du1;

		_K.push_back(T(cross_index0, cross_index1, Fu0du1));
		_K.push_back(T(cross_index1, cross_index0, Fu1du0));
	}
}

void ParallelContactSpring::solveV(vector<T>& _K, VectorXd& f, int nodes_size)
{
	double l = (node1->position - node0->position).norm();
	double delta_v = abs(node1->v - node0->v);

	if (delta_v >= d) return;

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d w = (node1->position - node0->position) / delta_v;

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d1 * d1.transpose();

	parallelContactEnergy = 0.5*Kc*L*(delta_v - d)*(delta_v - d);

	if (!node0->onBorder)
	{
		int cross_index0 = nodes_size * 3 + node0->cross_index * 2;
		double Fv0 = Kc * L*(delta_v - d);
		f(cross_index0 + 1) += Fv0;

		double Fv0dv0 = -Kc * L;
		_K.push_back(T(cross_index0 + 1, cross_index0 + 1, Fv0dv0));
	}

	if (!node1->onBorder)
	{
		int cross_index1 = nodes_size * 3 + node1->cross_index * 2;
		double Fv1 = -Kc * L*(delta_v - d);
		f(cross_index1 + 1) += Fv1;

		double Fv1dv1 = -Kc * L;
		_K.push_back(T(cross_index1 + 1, cross_index1 + 1, Fv1dv1));
	}

	if (!node0->onBorder && !node1->onBorder)
	{
		int cross_index0 = nodes_size * 3 + node0->cross_index * 2;
		int cross_index1 = nodes_size * 3 + node1->cross_index * 2;

		double Fv0dv1 = Kc * L;
		double Fv1dv0 = Fv0dv1;

		_K.push_back(T(cross_index0 + 1, cross_index1 + 1, Fv0dv1));
		_K.push_back(T(cross_index1 + 1, cross_index0 + 1, Fv1dv0));
	}
}