#include "bend.h"

using namespace std;
using namespace Eigen;

BendSpring::BendSpring(Node* n0, Node* n1, Node* n2, double _B, double _R, YarnType type)
{
	node0 = n0;
	
	if (type == Warp)
	{
		node1 = n1->u > n2->u ? n1 : n2;
		node2 = n2->u < n1->u ? n2 : n1;
	}
	else if (type == Weft)
	{
		node1 = n1->v > n2->v ? n1 : n2;
		node2 = n2->v < n1->v ? n2 : n1;
	}

	B = _B;
	R = _R;
	Kb = B * M_PI*R*R;

	bendEnergy = 0;

	springType = type;
}

BendSpring::~BendSpring() {}

void BendSpring::solve(vector<T>& _K, VectorXd& f)
{
	if (springType == Warp) solveU(_K, f);

	else if (springType == Weft) solveV(_K, f);

	return;
}


void BendSpring::solveU(vector<T>& _K, VectorXd& f)
{
	/*Offset to get crimp and then compute bending force*/
	Vector3d pos0 = node0->position - R * node0->normal;
	Vector3d pos1 = node1->position - R * node1->normal;
	Vector3d pos2 = node2->position - R * node2->normal;

	double l1 = (pos1 - pos0).norm();
	double l2 = (pos2 - pos0).norm();

	Vector3d d1 = pos1 - pos0; d1.normalize();
	Vector3d d2 = pos2 - pos0; d2.normalize();

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P1 = I - d1 * d1.transpose();
	Matrix3d P2 = I - d2 * d2.transpose();

	double theta = acos(-d1.dot(d2));
	
	double u0 = node0->u;
	double u1 = node1->u;
	double u2 = node2->u;

	bendEnergy = Kb * (theta*theta) / abs(u1 - u2);
	
	int index0 = node0->index * 5;
	int index1 = node1->index * 5;
	int index2 = node2->index * 5;

	//Compute and fill the force vector
	
	Vector3d Fq1 = (-2 * Kb*theta) / (l1*(u1 - u2)*sin(theta))*P1*d2;
	Vector3d Fq2 = (-2 * Kb*theta) / (l2*(u1 - u2)*sin(theta))*P2*d1;
	Vector3d Fq0 = -(Fq1 + Fq2);

	double Fu1 = Kb * theta*theta / ((u1 - u2)*(u1 - u2));
	double Fu2 = -Fu1;
	double Fu0 = 0;

	f.segment<3>(index0) += Fq0;
	f(index0 + 3) += Fu0;
	
	f.segment<3>(index1) += Fq1;
	f(index1 + 3) += Fu1;

	f.segment<3>(index2) += Fq2;
	f(index2 + 3) += Fu2;

	node0->compressForce += 0.5*node0->normal.dot(Fq0);
	node1->compressForce += 0.5*node1->normal.dot(Fq1);
	node2->compressForce += 0.5*node2->normal.dot(Fq2);
	
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
	
	MatrixXd f0q0;
	f0q0.resize(5, 5);
	f0q0.setZero();
	f0q0.block<3, 3>(0, 0) = Fq0dq0;
	fillBlock(_K, f0q0, index0, index0);

	MatrixXd f0q1;
	f0q1.resize(5, 5);
	f0q1.setZero();
	f0q1.block<3, 3>(0, 0) = Fq0dq1;
	f0q1.block<3, 1>(0, 3) = Fq0du1;
	fillBlock(_K, f0q1, index0, index1);

	MatrixXd f0q2;
	f0q2.resize(5, 5);
	f0q2.setZero();
	f0q2.block<3, 3>(0, 0) = Fq0dq2;
	f0q2.block<3, 1>(0, 3) = Fq0du2;
	fillBlock(_K, f0q2, index0, index2);

	MatrixXd f1q0;
	f1q0.resize(5, 5);
	f1q0.setZero();
	f1q0.block<3, 3>(0, 0) = Fq1dq0;
	f1q0.block<1, 3>(3, 0) = Fu1dq0;
	fillBlock(_K, f1q0, index1, index0);

	MatrixXd f1q1;
	f1q1.resize(5, 5);
	f1q1.setZero();
	f1q1.block<3, 3>(0, 0) = Fq1dq1;
	f1q1.block<3, 1>(0, 3) = Fq1du1;
	f1q1.block<1, 3>(3, 0) = Fu1dq1;
	f1q1(3, 3) = Fu1du1;
	fillBlock(_K, f1q1, index1, index1);

	MatrixXd f1q2;
	f1q2.resize(5, 5);
	f1q2.setZero();
	f1q2.block<3, 3>(0, 0) = Fq1dq2;
	f1q2.block<3, 1>(0, 3) = Fq1du2;
	f1q2.block<1, 3>(3, 0) = Fu1dq2;
	f1q2(3, 3) = Fu1du2;
	fillBlock(_K, f1q2, index1, index2);

	MatrixXd f2q0;
	f2q0.resize(5, 5);
	f2q0.setZero();
	f2q0.block<3, 3>(0, 0) = Fq2dq0;
	f2q0.block<1, 3>(3, 0) = Fu2dq0;
	fillBlock(_K, f2q0, index2, index0);

	MatrixXd f2q1;
	f2q1.resize(5, 5);
	f2q1.setZero();
	f2q1.block<3, 3>(0, 0) = Fq2dq1;
	f2q1.block<3, 1>(0, 3) = Fq2du1;
	f2q1.block<1, 3>(3, 0) = Fu2dq1;
	f2q1(3, 3) = Fu2du1;
	fillBlock(_K, f2q1, index2, index1);

	MatrixXd f2q2;
	f2q2.resize(5, 5);
	f2q2.setZero();
	f2q2.block<3, 3>(0, 0) = Fq2dq2;
	f2q2.block<3, 1>(0, 3) = Fq2du2;
	f2q2.block<1, 3>(3, 0) = Fu2dq2;
	f2q2(3, 3) = Fu2du2;
	fillBlock(_K, f2q2, index2, index2);

}

void BendSpring::solveV(vector<T>& _K, VectorXd& f)
{
	/*Offset to get crimp and then compute bending force*/
	Vector3d pos0 = node0->position + R * node0->normal;
	Vector3d pos1 = node1->position + R * node1->normal;
	Vector3d pos2 = node2->position + R * node2->normal;

	double l1 = (pos1 - pos0).norm();
	double l2 = (pos2 - pos0).norm();

	Vector3d d1 = pos1 - pos0; d1.normalize();
	Vector3d d2 = pos2 - pos0; d2.normalize();

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P1 = I - d1 * d1.transpose();
	Matrix3d P2 = I - d2 * d2.transpose();

	double theta = acos(-d1.dot(d2));

	double v0 = node0->v;
	double v1 = node1->v;
	double v2 = node2->v;

	bendEnergy = Kb * (theta*theta) / abs(v1 - v2);

	int index0 = node0->index * 5;
	int index1 = node1->index * 5;
	int index2 = node2->index * 5;

	//Compute and fill the force vector

	Vector3d Fq1 = (-2 * Kb*theta) / (l1*(v1 - v2)*sin(theta))*P1*d2;
	Vector3d Fq2 = (-2 * Kb*theta) / (l2*(v1 - v2)*sin(theta))*P2*d1;
	Vector3d Fq0 = -(Fq1 + Fq2);

	double Fv1 = Kb * theta*theta / ((v1 - v2)*(v1 - v2));
	double Fv2 = -Fv1;
	double Fv0 = 0;

	f.segment<3>(index0) += Fq0;
	f(index0 + 4) += Fv0;

	f.segment<3>(index1) += Fq1;
	f(index1 + 4) += Fv1;

	f.segment<3>(index2) += Fq2;
	f(index2 + 4) += Fv2;

	node0->compressForce -= 0.5*node0->normal.dot(Fq0);
	node1->compressForce -= 0.5*node1->normal.dot(Fq1);
	node2->compressForce -= 0.5*node2->normal.dot(Fq2);

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

	MatrixXd f0q0;
	f0q0.resize(5, 5);
	f0q0.setZero();
	f0q0.block<3, 3>(0, 0) = Fq0dq0;
	fillBlock(_K, f0q0, index0, index0);

	MatrixXd f0q1;
	f0q1.resize(5, 5);
	f0q1.setZero();
	f0q1.block<3, 3>(0, 0) = Fq0dq1;
	f0q1.block<3, 1>(0, 4) = Fq0dv1;
	fillBlock(_K, f0q1, index0, index1);

	MatrixXd f0q2;
	f0q2.resize(5, 5);
	f0q2.setZero();
	f0q2.block<3, 3>(0, 0) = Fq0dq2;
	f0q2.block<3, 1>(0, 4) = Fq0dv2;
	fillBlock(_K, f0q2, index0, index2);

	MatrixXd f1q0;
	f1q0.resize(5, 5);
	f1q0.setZero();
	f1q0.block<3, 3>(0, 0) = Fq1dq0;
	f1q0.block<1, 3>(4, 0) = Fv1dq0;
	fillBlock(_K, f1q0, index1, index0);

	MatrixXd f1q1;
	f1q1.resize(5, 5);
	f1q1.setZero();
	f1q1.block<3, 3>(0, 0) = Fq1dq1;
	f1q1.block<3, 1>(0, 4) = Fq1dv1;
	f1q1.block<1, 3>(4, 0) = Fv1dq1;
	f1q1(4, 4) = Fv1dv1;
	fillBlock(_K, f1q1, index1, index1);

	MatrixXd f1q2;
	f1q2.resize(5, 5);
	f1q2.setZero();
	f1q2.block<3, 3>(0, 0) = Fq1dq2;
	f1q2.block<3, 1>(0, 4) = Fq1dv2;
	f1q2.block<1, 3>(4, 0) = Fv1dq2;
	f1q2(4, 4) = Fv1dv2;
	fillBlock(_K, f1q2, index1, index2);

	MatrixXd f2q0;
	f2q0.resize(5, 5);
	f2q0.setZero();
	f2q0.block<3, 3>(0, 0) = Fq2dq0;
	f2q0.block<1, 3>(4, 0) = Fv2dq0;
	fillBlock(_K, f2q0, index2, index0);

	MatrixXd f2q1;
	f2q1.resize(5, 5);
	f2q1.setZero();
	f2q1.block<3, 3>(0, 0) = Fq2dq1;
	f2q1.block<3, 1>(0, 4) = Fq2dv1;
	f2q1.block<1, 3>(4, 0) = Fv2dq1;
	f2q1(4, 4) = Fv2dv1;
	fillBlock(_K, f2q1, index2, index1);

	MatrixXd f2q2;
	f2q2.resize(5, 5);
	f2q2.setZero();
	f2q2.block<3, 3>(0, 0) = Fq2dq2;
	f2q2.block<3, 1>(0, 4) = Fq2dv2;
	f2q2.block<1, 3>(4, 0) = Fv2dq2;
	f2q2(4, 4) = Fv2dv2;
	fillBlock(_K, f2q2, index2, index2);
}
