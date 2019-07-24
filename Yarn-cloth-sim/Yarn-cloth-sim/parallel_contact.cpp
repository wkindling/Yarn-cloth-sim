#include "parallel_contact.h"

using namespace std;
using namespace Eigen;

ParallelContactSpring::ParallelContactSpring(Node* n0, Node* n1, double _Kc, double _R, double _L, YarnType type)
{
	node0 = n0;
	node1 = n1;

	Kc = _Kc;
	L = _L;
	R = _R;

	d = 4 * R;
	parallelContactEnergy = 0;

	springType = type;
}

ParallelContactSpring::~ParallelContactSpring() {}

void ParallelContactSpring::solve()
{
	if (springType == Weft) solveV();

	else if (springType == Warp) solveU();
	
	return;
}

void ParallelContactSpring::solveU()
{
	double l = (node1->position - node0->position).norm();
	double delta_u = abs(node1->u - node0->u);

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d w = (node1->position - node0->position) / delta_u;

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d1 * d1.transpose();

	double V = 0.5*Kc*L*(delta_u - d)*(delta_u - d);

	//Compute and fill the force vector
	//This force will have no component on Lag position
	double Fu0 = Kc * L*(delta_u - d);
	double Fu1 = -Fu0;

	//Compute and fill the stiffness matrix
	//The local stiffness matrix should be 10*10
	double Fu0du0 = -Kc * L;
	double Fu1du1 = Fu0du0;
	double Fu0du1 = -Fu0du0;
	double Fu1du0 = Fu0du1;

	/* TODO: Fill the block just like EoL*/





}

void ParallelContactSpring::solveV()
{
	double l = (node1->position - node0->position).norm();
	double delta_v = abs(node1->v - node0->v);

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d w = (node1->position - node0->position) / delta_v;

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d1 * d1.transpose();

	double V = 0.5*Kc*L*(delta_v - d)*(delta_v - d);

	//Compute and fill the force vector
	//This force will have no component on Lag position
	double Fv0 = Kc * L*(delta_v - d);
	double Fv1 = -Fv0;

	//Compute and fill the stiffness matrix
	//The local stiffness matrix should be 10*10
	double Fv0dv0 = -Kc * L;
	double Fv1dv1 = Fv0dv0;
	double Fv0dv1 = -Fv0dv0;
	double Fv1dv0 = Fv0dv1;

	/* TODO: Fill the block just like EoL*/
}