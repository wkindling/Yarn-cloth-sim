#include "parallel_contact.h"

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

void ParallelContactSpring::solve(vector<T>& _K, VectorXd& f)
{
	if (springType == Weft) solveV(_K, f);

	else if (springType == Warp) solveU(_K, f);
	
	return;
}

void ParallelContactSpring::solveU(vector<T>& _K, VectorXd& f)
{
	double l = (node1->position - node0->position).norm();
	double delta_u = abs(node1->u - node0->u);

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d w = (node1->position - node0->position) / delta_u;

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d1 * d1.transpose();

	parallelContactEnergy = 0.5*Kc*L*(delta_u - d)*(delta_u - d);

	int index0 = node0->index * 5;
	int index1 = node1->index * 5;

	//Compute and fill the force vector
	//This force will have no component on Lag position
	double Fu0 = Kc * L*(delta_u - d);
	double Fu1 = -Fu0;

	f(index0 + 3) += Fu0;
	f(index1 + 3) += Fu1;

	//Compute and fill the stiffness matrix
	//The local stiffness matrix should be 10*10
	double Fu0du0 = -Kc * L;
	double Fu1du1 = Fu0du0;
	double Fu0du1 = -Fu0du0;
	double Fu1du0 = Fu0du1;

	_K.push_back(T(index0 + 3, index0 + 3, Fu0du0));
	_K.push_back(T(index0 + 3, index1 + 3, Fu0du1));
	_K.push_back(T(index1 + 3, index0 + 3, Fu1du0));
	_K.push_back(T(index1 + 3, index1 + 3, Fu1du1));
}

void ParallelContactSpring::solveV(vector<T>& _K, VectorXd& f)
{
	double l = (node1->position - node0->position).norm();
	double delta_v = abs(node1->v - node0->v);

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d w = (node1->position - node0->position) / delta_v;

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d1 * d1.transpose();

	parallelContactEnergy = 0.5*Kc*L*(delta_v - d)*(delta_v - d);

	int index0 = node0->index * 5;
	int index1 = node1->index * 5;

	//Compute and fill the force vector
	//This force will have no component on Lag position
	double Fv0 = Kc * L*(delta_v - d);
	double Fv1 = -Fv0;

	f(index0 + 4) += Fv0;
	f(index1 + 4) += Fv1;

	//Compute and fill the stiffness matrix
	//The local stiffness matrix should be 10*10
	double Fv0dv0 = -Kc * L;
	double Fv1dv1 = Fv0dv0;
	double Fv0dv1 = -Fv0dv0;
	double Fv1dv0 = Fv0dv1;
	
	_K.push_back(T(index0 + 4, index0 + 4, Fv0dv0));
	_K.push_back(T(index0 + 4, index1 + 4, Fv0dv1));
	_K.push_back(T(index1 + 4, index0 + 4, Fv1dv0));
	_K.push_back(T(index1 + 4, index1 + 4, Fv1dv1));
}