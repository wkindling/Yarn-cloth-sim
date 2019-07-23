#include "bend.h"

using namespace std;
using namespace Eigen;

BendSpring::BendSpring(Node* n0, Node* n1, Node* n2, double _B, double _R, YarnType type)
{
	node0 = n0;
	node1 = n1;
	node2 = n2;

	B = _B;
	R = _R;
	Kb = B * M_PI*R*R;

	bendEnergy = 0;

	springType = type;
}

BendSpring::~BendSpring() {}

void BendSpring::solve()
{
	if (springType == Warp) solveU();

	else if (springType == Weft) solveV();

	return;
}


void BendSpring::solveU()
{
	double l1 = (node1->position - node0->position).norm();
	double l2 = (node2->position - node0->position).norm();

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d d2 = node2->position - node0->position; d2.normalize();

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P1 = I - d1 * d1.transpose();
	Matrix3d P2 = I - d2 * d2.transpose();

	double theta = acos(-d1.dot(d2));
	
	double u0 = node0->u;
	double u1 = node1->u;
	double u2 = node2->u;
	
	//Compute and fill the force vector
	
	Vector3d Fq1 = (-2 * Kb*theta) / (l1*(u1 - u2)*sin(theta))*P1*d2;
	Vector3d Fq2 = (-2 * Kb*theta) / (l2*(u1 - u2)*sin(theta))*P2*d1;
	Vector3d Fq0 = -(Fq1 + Fq2);

	double Fu1 = Kb * theta*theta / ((u1 - u2)*(u1 - u2));
	double Fu2 = -Fu1;
	double Fu0 = 0;

	//Compute and fill the stiffness matrix
	//The local stiffness matrix should be 15*15
	Matrix3d Fq1dq1 = 2 * Kb / (l1*l1*(u1 - u2)*sin(theta))*(theta*(P1*d2*d1.transpose() + cos(theta) / sin(theta) / sin(theta)*P1*d2*d2.transpose()*P1 + cos(theta)*P1 + d1 * d2.transpose()*P1) - 1 / (sin(theta))*P1*d2*d2.transpose()*P1);
	Matrix3d Fq1dq2 = -(2 * Kb) / (l2*l1*(u1 - u2)*sin(theta))*(theta*(P1 - cos(theta) / sin(theta) / sin(theta)*P1*d2*d1.transpose()) + 1 / (sin(theta))*P1*d2*d1.transpose())*P2;
	Matrix3d Fq2dq1 = -(2 * Kb) / (l1*l2*(u1 - u2)*sin(theta))*(theta*(P2 - cos(theta) / sin(theta) / sin(theta)*P2*d1*d2.transpose()) + 1 / (sin(theta))*P2*d1*d2.transpose())*P1;
	Matrix3d Fq2dq2 = 2 * Kb / (l2*l2*(u1 - u2)*sin(theta))*(theta*(P2*d1*d2.transpose() + cos(theta) / sin(theta) / sin(theta)*P2*d1*d1.transpose()*P2 + cos(theta)*P2 + d2 * d1.transpose()*P2) - 1 / (sin(theta))*P2*d1*d1.transpose()*P2);

	Matrix3d Fq1dq0 = -(Fq1dq1 + Fq1dq2);
	Matrix3d Fq2dq0 = -(Fq2dq1 + Fq2dq2);
	Matrix3d Fq0dq1 = -(Fq1dq1 + Fq2dq1);
	Matrix3d Fq0dq2 = -(Fq1dq2 + Fq2dq2);
	Matrix3d Fq0dq0 = -(Fq1dq0 + Fq2dq0);

	double Fu1du1 = -2 * Kb*theta*theta / pow((u1 - u2), 3);
	double Fu2du2 = Fu1du1;
	double Fu1du2 = -Fu1du1;
	double Fu2du1 = Fu1du2;

	Vector3d Fq1du1 = 2 * Kb*theta / (l1*(u1 - u2)*(u1 - u2)*sin(theta))*P1*d2;
	Vector3d Fq1du2 = -Fq1du1;
	Vector3d Fq2du1 = 2 * Kb*theta / (l2*(u1 - u2)*(u1 - u2)*sin(theta))*P2*d1;
	Vector3d Fq2du2 = -Fq2du1;

	Vector3d Fq0du1 = -(Fq1du1 + Fq2du1);
	Vector3d Fq0du2 = -Fq0du1;

	Vector3d Fu1dq1 = 2 * Kb*theta / (l1*(u1 - u2)*(u1 - u2)*sin(theta))*d2.transpose()*P1;
	Vector3d Fu2dq1 = -Fu1dq1;

	Vector3d Fu1dq2 = 2 * Kb*theta / (l2*(u1 - u2)*(u1 - u2)*sin(theta))*d1.transpose()*P2;
	Vector3d Fu2dq2 = -Fu1dq2;

	Vector3d Fu1dq0 = -(Fu1dq1 + Fu1dq2);
	Vector3d Fu2dq0 = -Fu1dq0;
	

	/* TODO: Fill the block just like EoL*/





}

void BendSpring::solveV()
{
	double l1 = (node1->position - node0->position).norm();
	double l2 = (node2->position - node0->position).norm();

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d d2 = node2->position - node0->position; d2.normalize();

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P1 = I - d1 * d1.transpose();
	Matrix3d P2 = I - d2 * d2.transpose();

	double theta = acos(-d1.dot(d2));

	double v0 = node0->v;
	double v1 = node1->v;
	double v2 = node2->v;

	//Compute and fill the force vector

	Vector3d Fq1 = (-2 * Kb*theta) / (l1*(v1 - v2)*sin(theta))*P1*d2;
	Vector3d Fq2 = (-2 * Kb*theta) / (l2*(v1 - v2)*sin(theta))*P2*d1;
	Vector3d Fq0 = -(Fq1 + Fq2);

	double Fv1 = Kb * theta*theta / ((v1 - v2)*(v1 - v2));
	double Fv2 = -Fv1;
	double Fv0 = 0;

	//Compute and fill the stiffness matrix
	//The local stiffness matrix should be 15*15
	Matrix3d Fq1dq1 = 2 * Kb / (l1*l1*(v1 - v2)*sin(theta))*(theta*(P1*d2*d1.transpose() + cos(theta) / sin(theta) / sin(theta)*P1*d2*d2.transpose()*P1 + cos(theta)*P1 + d1 * d2.transpose()*P1) - 1 / (sin(theta))*P1*d2*d2.transpose()*P1);
	Matrix3d Fq1dq2 = -(2 * Kb) / (l2*l1*(v1 - v2)*sin(theta))*(theta*(P1 - cos(theta) / sin(theta) / sin(theta)*P1*d2*d1.transpose()) + 1 / (sin(theta))*P1*d2*d1.transpose())*P2;
	Matrix3d Fq2dq1 = -(2 * Kb) / (l1*l2*(v1 - v2)*sin(theta))*(theta*(P2 - cos(theta) / sin(theta) / sin(theta)*P2*d1*d2.transpose()) + 1 / (sin(theta))*P2*d1*d2.transpose())*P1;
	Matrix3d Fq2dq2 = 2 * Kb / (l2*l2*(v1 - v2)*sin(theta))*(theta*(P2*d1*d2.transpose() + cos(theta) / sin(theta) / sin(theta)*P2*d1*d1.transpose()*P2 + cos(theta)*P2 + d2 * d1.transpose()*P2) - 1 / (sin(theta))*P2*d1*d1.transpose()*P2);

	Matrix3d Fq1dq0 = -(Fq1dq1 + Fq1dq2);
	Matrix3d Fq2dq0 = -(Fq2dq1 + Fq2dq2);
	Matrix3d Fq0dq1 = -(Fq1dq1 + Fq2dq1);
	Matrix3d Fq0dq2 = -(Fq1dq2 + Fq2dq2);
	Matrix3d Fq0dq0 = -(Fq1dq0 + Fq2dq0);

	double Fv1dv1 = -2 * Kb*theta*theta / pow((v1 - v2), 3);
	double Fv2dv2 = Fv1dv1;
	double Fv1dv2 = -Fv1dv1;
	double Fv2dv1 = Fv1dv2;

	Vector3d Fq1dv1 = 2 * Kb*theta / (l1*(v1 - v2)*(v1 - v2)*sin(theta))*P1*d2;
	Vector3d Fq1dv2 = -Fq1dv1;
	Vector3d Fq2dv1 = 2 * Kb*theta / (l2*(v1 - v2)*(v1 - v2)*sin(theta))*P2*d1;
	Vector3d Fq2dv2 = -Fq2dv1;

	Vector3d Fq0dv1 = -(Fq1dv1 + Fq2dv1);
	Vector3d Fq0dv2 = -Fq0dv1;

	Vector3d Fv1dq1 = 2 * Kb*theta / (l1*(v1 - v2)*(v1 - v2)*sin(theta))*d2.transpose()*P1;
	Vector3d Fv2dq1 = -Fv1dq1;

	Vector3d Fv1dq2 = 2 * Kb*theta / (l2*(v1 - v2)*(v1 - v2)*sin(theta))*d1.transpose()*P2;
	Vector3d Fv2dq2 = -Fv1dq2;

	Vector3d Fv1dq0 = -(Fv1dq1 + Fv1dq2);
	Vector3d Fv2dq0 = -Fv1dq0;


	/* TODO: Fill the block just like EoL*/




}
