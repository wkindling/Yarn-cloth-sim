#include "bend.h"
#include <iostream>

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

void BendSpring::solve(vector<T>& _K, VectorXd& f, int nodes_size, double h)
{
	if (springType == Warp) solveU(_K, f, nodes_size, h);

	else if (springType == Weft) solveV(_K, f, nodes_size, h);

	return;
}


void BendSpring::solveU(vector<T>& _K, VectorXd& f, int nodes_size, double h)
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
	
	int node_index0 = node0->node_index * 3;
	int node_index1 = node1->node_index * 3;
	int node_index2 = node2->node_index * 3;

	int cross_index0 = nodes_size * 3 + node0->cross_index * 2;
	int cross_index1 = nodes_size * 3 + node1->cross_index * 2;
	int cross_index2 = nodes_size * 3 + node2->cross_index * 2;

	//Compute and fill the force vector
	Vector3d Fq1 = (-2 * Kb*theta) / (l1*(u1 - u2)*sin(theta))*P1*d2;
	Vector3d Fq2 = (-2 * Kb*theta) / (l2*(u1 - u2)*sin(theta))*P2*d1;
	Vector3d Fq0 = -(Fq1 + Fq2);

	f.segment<3>(node_index0) += Fq0;
	f.segment<3>(node_index1) += Fq1;
	f.segment<3>(node_index2) += Fq2;

	node0->compressForce += 0.5*node0->normal.dot(Fq0);
	node1->compressForce += 0.5*node1->normal.dot(Fq1);
	node2->compressForce += 0.5*node2->normal.dot(Fq2);

	if (!node1->onBorder)
	{
		double Fu1 = Kb * theta*theta / ((u1 - u2)*(u1 - u2));
		f(cross_index1) += Fu1;
	}
	if (!node2->onBorder)
	{
		double Fu2 = -Kb * theta*theta / ((u1 - u2)*(u1 - u2));
		f(cross_index2) += Fu2;
	}

	//Compute and fill the stiffness matrix
	/* Lagrange Part */
	MatrixXd Fq1dq1 = 2 * Kb / (l1*l1*(u1 - u2)*sin(theta))*(theta*(P1*d2*d1.transpose() + cos(theta) / sin(theta) / sin(theta)*P1*d2*d2.transpose()*P1 + cos(theta)*P1 + d1 * d2.transpose()*P1) - 1 / (sin(theta))*P1*d2*d2.transpose()*P1);
	MatrixXd Fq1dq2 = -(2 * Kb) / (l2*l1*(u1 - u2)*sin(theta))*(theta*(P1 - cos(theta) / sin(theta) / sin(theta)*P1*d2*d1.transpose()) + 1 / (sin(theta))*P1*d2*d1.transpose())*P2;
	MatrixXd Fq2dq1 = -(2 * Kb) / (l1*l2*(u1 - u2)*sin(theta))*(theta*(P2 - cos(theta) / sin(theta) / sin(theta)*P2*d1*d2.transpose()) + 1 / (sin(theta))*P2*d1*d2.transpose())*P1;
	MatrixXd Fq2dq2 = 2 * Kb / (l2*l2*(u1 - u2)*sin(theta))*(theta*(P2*d1*d2.transpose() + cos(theta) / sin(theta) / sin(theta)*P2*d1*d1.transpose()*P2 + cos(theta)*P2 + d2 * d1.transpose()*P2) - 1 / (sin(theta))*P2*d1*d1.transpose()*P2);

	MatrixXd Fq1dq0 = -(Fq1dq1 + Fq1dq2);
	MatrixXd Fq2dq0 = -(Fq2dq1 + Fq2dq2);
	MatrixXd Fq0dq1 = -(Fq1dq1 + Fq2dq1);
	MatrixXd Fq0dq2 = -(Fq1dq2 + Fq2dq2);
	MatrixXd Fq0dq0 = -(Fq1dq0 + Fq2dq0);

	fillGlobalStiffness(_K, Fq1dq1, node_index1, node_index1, h);
	fillGlobalStiffness(_K, Fq1dq2, node_index1, node_index2, h);
	fillGlobalStiffness(_K, Fq1dq0, node_index1, node_index0, h);
	fillGlobalStiffness(_K, Fq2dq1, node_index2, node_index1, h);
	fillGlobalStiffness(_K, Fq2dq2, node_index2, node_index2, h);
	fillGlobalStiffness(_K, Fq2dq0, node_index2, node_index0, h);	
	fillGlobalStiffness(_K, Fq0dq1, node_index0, node_index1, h);
	fillGlobalStiffness(_K, Fq0dq2, node_index0, node_index2, h);
	fillGlobalStiffness(_K, Fq0dq0, node_index0, node_index0, h);

	if (!node0->onBorder)
	{
		MatrixXd fri_u0dq0 = 0.5*node0->normal.transpose()*Fq0dq0;
		MatrixXd fri_u0dq1 = 0.5*node0->normal.transpose()*Fq0dq1;
		MatrixXd fri_u0dq2 = 0.5*node0->normal.transpose()*Fq0dq2;

		fillGlobalStiffness(node0->u_friction, fri_u0dq0, cross_index0, node_index0, h);
		fillGlobalStiffness(node0->u_friction, fri_u0dq1, cross_index0, node_index1, h);
		fillGlobalStiffness(node0->u_friction, fri_u0dq2, cross_index0, node_index2, h);
	}

	if (!node1->onBorder)
	{
		MatrixXd fri_u1dq0 = 0.5*node1->normal.transpose()*Fq1dq0;
		MatrixXd fri_u1dq1 = 0.5*node1->normal.transpose()*Fq1dq1;
		MatrixXd fri_u1dq2 = 0.5*node1->normal.transpose()*Fq1dq2;

		fillGlobalStiffness(node1->u_friction, fri_u1dq0, cross_index1, node_index0, h);
		fillGlobalStiffness(node1->u_friction, fri_u1dq1, cross_index1, node_index1, h);
		fillGlobalStiffness(node1->u_friction, fri_u1dq2, cross_index1, node_index2, h);
	}

	if (!node2->onBorder)
	{
		MatrixXd fri_u2dq0 = 0.5*node2->normal.transpose()*Fq2dq0;
		MatrixXd fri_u2dq1 = 0.5*node2->normal.transpose()*Fq2dq1;
		MatrixXd fri_u2dq2 = 0.5*node2->normal.transpose()*Fq2dq2;

		fillGlobalStiffness(node2->u_friction, fri_u2dq0, cross_index2, node_index0, h);
		fillGlobalStiffness(node2->u_friction, fri_u2dq1, cross_index2, node_index1, h);
		fillGlobalStiffness(node2->u_friction, fri_u2dq2, cross_index2, node_index2, h);
	}

	/* Euler Part */
	if (!node1->onBorder)
	{
		double Fu1du1 = -2 * Kb*theta*theta / pow((u1 - u2), 3);
		_K.push_back(T(cross_index1, cross_index1, -h * h*Fu1du1));

		MatrixXd Fq1du1 = 2 * Kb*theta / (l1*(u1 - u2)*(u1 - u2)*sin(theta))*P1*d2;
		MatrixXd Fq2du1 = 2 * Kb*theta / (l2*(u1 - u2)*(u1 - u2)*sin(theta))*P2*d1;
		MatrixXd Fq0du1 = -(Fq1du1 + Fq2du1);
		fillGlobalStiffness(_K, Fq1du1, node_index1, cross_index1, h);
		fillGlobalStiffness(_K, Fq2du1, node_index2, cross_index1, h);
		fillGlobalStiffness(_K, Fq0du1, node_index0, cross_index1, h);

		if (!node0->onBorder)
		{
			MatrixXd fri_u0du1 = 0.5*node0->normal.transpose()*Fq0du1;
			fillGlobalStiffness(node0->u_friction, fri_u0du1, cross_index0, cross_index1, h);
		}

		if (!node1->onBorder)
		{
			MatrixXd fri_u1du1 = 0.5*node1->normal.transpose()*Fq1du1;
			fillGlobalStiffness(node1->u_friction, fri_u1du1, cross_index1, cross_index1, h);
		}

		if (!node2->onBorder)
		{
			MatrixXd fri_u2du1 = 0.5*node2->normal.transpose()*Fq2du1;
			fillGlobalStiffness(node2->u_friction, fri_u2du1, cross_index2, cross_index1, h);
		}

		MatrixXd Fu1dq1 = 2 * Kb*theta / (l1*(u1 - u2)*(u1 - u2)*sin(theta))*d2.transpose()*P1;
		MatrixXd Fu1dq2 = 2 * Kb*theta / (l2*(u1 - u2)*(u1 - u2)*sin(theta))*d1.transpose()*P2;
		MatrixXd Fu1dq0 = -(Fu1dq1 + Fu1dq2);
		fillGlobalStiffness(_K, Fu1dq1, cross_index1, node_index1, h);
		fillGlobalStiffness(_K, Fu1dq2, cross_index1, node_index2, h);
		fillGlobalStiffness(_K, Fu1dq0, cross_index1, node_index0, h);
	}

	if (!node2->onBorder)
	{
		double Fu2du2 = -2 * Kb*theta*theta / pow((u1 - u2), 3);
		_K.push_back(T(cross_index2, cross_index2, -h*h*Fu2du2));

		MatrixXd Fq1du2 = -2 * Kb*theta / (l1*(u1 - u2)*(u1 - u2)*sin(theta))*P1*d2;
		MatrixXd Fq2du2 = -2 * Kb*theta / (l2*(u1 - u2)*(u1 - u2)*sin(theta))*P2*d1;
		MatrixXd Fq0du2 = -(Fq1du2 + Fq2du2);
		fillGlobalStiffness(_K, Fq1du2, node_index1, cross_index2, h);
		fillGlobalStiffness(_K, Fq2du2, node_index2, cross_index2, h);
		fillGlobalStiffness(_K, Fq0du2, node_index0, cross_index2, h);
		
		if (!node0->onBorder)
		{
			MatrixXd fri_u0du2 = 0.5*node0->normal.transpose()*Fq0du2;
			fillGlobalStiffness(node0->u_friction, fri_u0du2, cross_index0, cross_index2, h);
		}

		if (!node1->onBorder)
		{
			MatrixXd fri_u1du2 = 0.5*node1->normal.transpose()*Fq1du2;
			fillGlobalStiffness(node1->u_friction, fri_u1du2, cross_index1, cross_index2, h);
		}

		if (!node2->onBorder)
		{
			MatrixXd fri_u2du2 = 0.5*node2->normal.transpose()*Fq2du2;
			fillGlobalStiffness(node2->u_friction, fri_u2du2, cross_index2, cross_index2, h);
		}

		MatrixXd Fu2dq1 = -2 * Kb*theta / (l1*(u1 - u2)*(u1 - u2)*sin(theta))*d2.transpose()*P1;
		MatrixXd Fu2dq2 = -2 * Kb*theta / (l2*(u1 - u2)*(u1 - u2)*sin(theta))*d1.transpose()*P2;
		MatrixXd Fu2dq0 = -(Fu2dq1 + Fu2dq2);
		fillGlobalStiffness(_K, Fu2dq1, cross_index2, node_index1, h);
		fillGlobalStiffness(_K, Fu2dq2, cross_index2, node_index2, h);
		fillGlobalStiffness(_K, Fu2dq0, cross_index2, node_index0, h);
	}

	if (!node1->onBorder && !node2->onBorder)
	{
		double Fu1du2 = 2 * Kb*theta*theta / pow((u1 - u2), 3);
		double Fu2du1 = Fu1du2;

		_K.push_back(T(cross_index1, cross_index2, -h*h*Fu1du2));
		_K.push_back(T(cross_index2, cross_index1, -h*h*Fu2du1));
	}
}

void BendSpring::solveV(vector<T>& _K, VectorXd& f, int nodes_size, double h)
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

	int node_index0 = node0->node_index * 3;
	int node_index1 = node1->node_index * 3;
	int node_index2 = node2->node_index * 3;

	int cross_index0 = nodes_size * 3 + node0->cross_index * 2;
	int cross_index1 = nodes_size * 3 + node1->cross_index * 2;
	int cross_index2 = nodes_size * 3 + node2->cross_index * 2;

	//Compute and fill the force vector

	Vector3d Fq1 = (-2 * Kb*theta) / (l1*(v1 - v2)*sin(theta))*P1*d2;
	Vector3d Fq2 = (-2 * Kb*theta) / (l2*(v1 - v2)*sin(theta))*P2*d1;
	Vector3d Fq0 = -(Fq1 + Fq2);

	f.segment<3>(node_index0) += Fq0;
	f.segment<3>(node_index1) += Fq1;
	f.segment<3>(node_index2) += Fq2;

	node0->compressForce -= 0.5*node0->normal.dot(Fq0);
	node1->compressForce -= 0.5*node1->normal.dot(Fq1);
	node2->compressForce -= 0.5*node2->normal.dot(Fq2);

	if (!node1->onBorder)
	{
		double Fv1 = Kb * theta*theta / ((v1 - v2)*(v1 - v2));
		f(cross_index1 + 1) += Fv1;
	}

	if (!node2->onBorder)
	{
		double Fv2 = -Kb * theta*theta / ((v1 - v2)*(v1 - v2));
		f(cross_index2 + 1) += Fv2;
	}

	//Compute and fill the stiffness matrix
	/* Lagrange Part */
	MatrixXd Fq1dq1 = 2 * Kb / (l1*l1*(v1 - v2)*sin(theta))*(theta*(P1*d2*d1.transpose() + cos(theta) / sin(theta) / sin(theta)*P1*d2*d2.transpose()*P1 + cos(theta)*P1 + d1 * d2.transpose()*P1) - 1 / (sin(theta))*P1*d2*d2.transpose()*P1);
	MatrixXd Fq1dq2 = -(2 * Kb) / (l2*l1*(v1 - v2)*sin(theta))*(theta*(P1 - cos(theta) / sin(theta) / sin(theta)*P1*d2*d1.transpose()) + 1 / (sin(theta))*P1*d2*d1.transpose())*P2;
	MatrixXd Fq2dq1 = -(2 * Kb) / (l1*l2*(v1 - v2)*sin(theta))*(theta*(P2 - cos(theta) / sin(theta) / sin(theta)*P2*d1*d2.transpose()) + 1 / (sin(theta))*P2*d1*d2.transpose())*P1;
	MatrixXd Fq2dq2 = 2 * Kb / (l2*l2*(v1 - v2)*sin(theta))*(theta*(P2*d1*d2.transpose() + cos(theta) / sin(theta) / sin(theta)*P2*d1*d1.transpose()*P2 + cos(theta)*P2 + d2 * d1.transpose()*P2) - 1 / (sin(theta))*P2*d1*d1.transpose()*P2);

	MatrixXd Fq1dq0 = -(Fq1dq1 + Fq1dq2);
	MatrixXd Fq2dq0 = -(Fq2dq1 + Fq2dq2);
	MatrixXd Fq0dq1 = -(Fq1dq1 + Fq2dq1);
	MatrixXd Fq0dq2 = -(Fq1dq2 + Fq2dq2);
	MatrixXd Fq0dq0 = -(Fq1dq0 + Fq2dq0);

	fillGlobalStiffness(_K, Fq1dq1, node_index1, node_index1, h);
	fillGlobalStiffness(_K, Fq1dq2, node_index1, node_index2, h);
	fillGlobalStiffness(_K, Fq1dq0, node_index1, node_index0, h);
	fillGlobalStiffness(_K, Fq2dq1, node_index2, node_index1, h);
	fillGlobalStiffness(_K, Fq2dq2, node_index2, node_index2, h);
	fillGlobalStiffness(_K, Fq2dq0, node_index2, node_index0, h);
	fillGlobalStiffness(_K, Fq0dq1, node_index0, node_index1, h);
	fillGlobalStiffness(_K, Fq0dq2, node_index0, node_index2, h);
	fillGlobalStiffness(_K, Fq0dq0, node_index0, node_index0, h);

	if (!node0->onBorder)
	{
		MatrixXd fri_v0dq0 = -0.5*node0->normal.transpose()*Fq0dq0;
		MatrixXd fri_v0dq1 = -0.5*node0->normal.transpose()*Fq0dq1;
		MatrixXd fri_v0dq2 = -0.5*node0->normal.transpose()*Fq0dq2;
		
		fillGlobalStiffness(node0->v_friction, fri_v0dq0, cross_index0 + 1, node_index0, h);
		fillGlobalStiffness(node0->v_friction, fri_v0dq1, cross_index0 + 1, node_index1, h);
		fillGlobalStiffness(node0->v_friction, fri_v0dq2, cross_index0 + 1, node_index2, h);
	}

	if (!node1->onBorder)
	{
		MatrixXd fri_v1dq0 = -0.5*node1->normal.transpose()*Fq1dq0;
		MatrixXd fri_v1dq1 = -0.5*node1->normal.transpose()*Fq1dq1;
		MatrixXd fri_v1dq2 = -0.5*node1->normal.transpose()*Fq1dq2;

		fillGlobalStiffness(node1->v_friction, fri_v1dq0, cross_index1 + 1, node_index0, h);
		fillGlobalStiffness(node1->v_friction, fri_v1dq1, cross_index1 + 1, node_index1, h);
		fillGlobalStiffness(node1->v_friction, fri_v1dq2, cross_index1 + 1, node_index2, h);
	}

	if (!node2->onBorder)
	{ 
		MatrixXd fri_v2dq0 = -0.5*node2->normal.transpose()*Fq2dq0;
		MatrixXd fri_v2dq1 = -0.5*node2->normal.transpose()*Fq2dq1;
		MatrixXd fri_v2dq2 = -0.5*node2->normal.transpose()*Fq2dq2;

		fillGlobalStiffness(node2->v_friction, fri_v2dq0, cross_index2 + 1, node_index0, h);
		fillGlobalStiffness(node2->v_friction, fri_v2dq1, cross_index2 + 1, node_index1, h);
		fillGlobalStiffness(node2->v_friction, fri_v2dq2, cross_index2 + 1, node_index2, h);
	}

	/* Euler Part */
	if (!node1->onBorder)
	{
		double Fv1dv1 = -2 * Kb*theta*theta / pow((v1 - v2), 3);
		_K.push_back(T(cross_index1 + 1, cross_index1 + 1, -h*h*Fv1dv1));

		MatrixXd Fq1dv1 = 2 * Kb*theta / (l1*(v1 - v2)*(v1 - v2)*sin(theta))*P1*d2;
		MatrixXd Fq2dv1 = 2 * Kb*theta / (l2*(v1 - v2)*(v1 - v2)*sin(theta))*P2*d1;
		MatrixXd Fq0dv1 = -(Fq1dv1 + Fq2dv1);
		fillGlobalStiffness(_K, Fq1dv1, node_index1, cross_index1 + 1, h);
		fillGlobalStiffness(_K, Fq2dv1, node_index2, cross_index1 + 1, h);
		fillGlobalStiffness(_K, Fq0dv1, node_index0, cross_index1 + 1, h);

		if (!node0->onBorder)
		{
			MatrixXd fri_v0dv1 = -0.5*node0->normal.transpose()*Fq0dv1;
			fillGlobalStiffness(node0->v_friction, fri_v0dv1, cross_index0 + 1, cross_index1 + 1, h);
		}

		if (!node1->onBorder)
		{
			MatrixXd fri_v1dv1 = -0.5*node1->normal.transpose()*Fq1dv1;
			fillGlobalStiffness(node1->v_friction, fri_v1dv1, cross_index1 + 1, cross_index1 + 1, h);
		}

		if (!node2->onBorder)
		{
			MatrixXd fri_v2dv1 = -0.5*node2->normal.transpose()*Fq2dv1;
			fillGlobalStiffness(node2->v_friction, fri_v2dv1, cross_index2 + 1, cross_index1 + 1, h);
		}

		MatrixXd Fv1dq1 = 2 * Kb*theta / (l1*(v1 - v2)*(v1 - v2)*sin(theta))*d2.transpose()*P1;
		MatrixXd Fv1dq2 = 2 * Kb*theta / (l2*(v1 - v2)*(v1 - v2)*sin(theta))*d1.transpose()*P2;
		MatrixXd Fv1dq0 = -(Fv1dq1 + Fv1dq2);
		fillGlobalStiffness(_K, Fv1dq1, cross_index1 + 1, node_index1, h);
		fillGlobalStiffness(_K, Fv1dq2, cross_index1 + 1, node_index2, h);
		fillGlobalStiffness(_K, Fv1dq0, cross_index1 + 1, node_index0, h);
	}

	if (!node2->onBorder)
	{
		double Fv2dv2 = -2 * Kb*theta*theta / pow((v1 - v2), 3);
		_K.push_back(T(cross_index2, cross_index2, -h*h*Fv2dv2));

		MatrixXd Fq1dv2 = -2 * Kb*theta / (l1*(v1 - v2)*(v1 - v2)*sin(theta))*P1*d2;
		MatrixXd Fq2dv2 = -2 * Kb*theta / (l2*(v1 - v2)*(v1 - v2)*sin(theta))*P2*d1;
		MatrixXd Fq0dv2 = -(Fq1dv2 + Fq2dv2);
		fillGlobalStiffness(_K, Fq1dv2, node_index1, cross_index2 + 1, h);
		fillGlobalStiffness(_K, Fq2dv2, node_index2, cross_index2 + 1, h);
		fillGlobalStiffness(_K, Fq0dv2, node_index0, cross_index2 + 1, h);

		if (!node0->onBorder)
		{
			MatrixXd fri_v0dv2 = -0.5*node0->normal.transpose()*Fq0dv2;
			fillGlobalStiffness(node0->v_friction, fri_v0dv2, cross_index0 + 1, cross_index2 + 1, h);
		}

		if (!node1->onBorder)
		{
			MatrixXd fri_v1dv2 = -0.5*node1->normal.transpose()*Fq1dv2;
			fillGlobalStiffness(node1->v_friction, fri_v1dv2, cross_index1 + 1, cross_index2 + 1, h);
		}

		if (!node2->onBorder)
		{
			MatrixXd fri_v2dv2 = -0.5*node2->normal.transpose()*Fq2dv2;
			fillGlobalStiffness(node2->v_friction, fri_v2dv2, cross_index2 + 1, cross_index2 + 1, h);
		}

		MatrixXd Fv2dq1 = -2 * Kb*theta / (l1*(v1 - v2)*(v1 - v2)*sin(theta))*d2.transpose()*P1;
		MatrixXd Fv2dq2 = -2 * Kb*theta / (l2*(v1 - v2)*(v1 - v2)*sin(theta))*d1.transpose()*P2;
		MatrixXd Fv2dq0 = -(Fv2dq1 + Fv2dq2);
		fillGlobalStiffness(_K, Fv2dq1, cross_index2 + 1, node_index1, h);
		fillGlobalStiffness(_K, Fv2dq2, cross_index2 + 1, node_index2, h);
		fillGlobalStiffness(_K, Fv2dq0, cross_index2 + 1, node_index0, h);
	}

	if (!node1->onBorder && !node2->onBorder)
	{
		double Fv1dv2 = 2 * Kb*theta*theta / pow((v1 - v2), 3);
		double Fv2dv1 = Fv1dv2;

		_K.push_back(T(cross_index1 + 1, cross_index2 + 1, -h*h*Fv1dv2));
		_K.push_back(T(cross_index2 + 1, cross_index1 + 1, -h*h*Fv2dv1));
	}

}
