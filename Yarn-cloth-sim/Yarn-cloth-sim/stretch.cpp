#include "stretch.h"
#include <freeglut.h>
#include <iostream>

using namespace std;
using namespace Eigen;

StretchSpring::StretchSpring(Node* n0, Node* n1, double _Y, double _R, YarnType type)
{
	//To ensure u1-u0>0 or v1-v0>0
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

	Y = _Y;
	R = _R;

	Ks = Y * M_PI*R*R;
	stretchEnergy = 0;
	springType = type;
}

StretchSpring::~StretchSpring() {}

void StretchSpring::solve(vector<T>& _K, VectorXd& f, int nodes_size, double h)
{
	if (springType == Weft) solveV(_K, f, nodes_size, h);
	
	else if (springType == Warp) solveU(_K, f, nodes_size, h);

	return;
}

void StretchSpring::solveU(vector<T>& _K, VectorXd& f, int nodes_size, double h)
{
	double l = (node1->position - node0->position).norm();
	double delta_u = abs(node1->u - node0->u);

	Vector3d w = (node1->position - node0->position) / delta_u;
	Vector3d d = node1->position - node0->position; d.normalize();
	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d * d.transpose();

	stretchEnergy = 0.5*Ks*delta_u*(w.norm() - 1)*(w.norm() - 1);

	int node_index0 = node0->node_index * 3;
	int node_index1 = node1->node_index * 3;

	int cross_index0 = nodes_size * 3 + node0->cross_index * 2;
	int cross_index1 = nodes_size * 3 + node1->cross_index * 2;

	/* Compute and fill the force vector */
	Vector3d Fq1 = -Ks * (w.norm() - 1)*d;
	Vector3d Fq0 = -Fq1;
	
	f.segment<3>(node_index0) += Fq0;
	f.segment<3>(node_index1) += Fq1;

	node0->compressForce += 0.5*node0->normal.dot(Fq0);
	node1->compressForce += 0.5*node1->normal.dot(Fq1);

	if (!node0->onBorder)
	{
		double Fu0 = -0.5*Ks*(w.norm()*w.norm() - 1);
		f(cross_index0) += Fu0;
	}

	if (!node1->onBorder)
	{
		double Fu1 = 0.5*Ks*(w.norm()*w.norm() - 1);
		f(cross_index1) += Fu1;
	}

	/* Compute and fill the stiffness matrix */
	/* Lagrange Part */
	MatrixXd Fq1dq1 = Ks / l * P - Ks / delta_u * I;
	MatrixXd Fq0dq0 = Fq1dq1;
	MatrixXd Fq1dq0 = -Fq1dq1;
	MatrixXd Fq0dq1 = Fq1dq0;

	fillGlobalStiffness(_K, Fq1dq1, node_index1, node_index1, h);
	fillGlobalStiffness(_K, Fq1dq0, node_index1, node_index0, h);
	fillGlobalStiffness(_K, Fq0dq0, node_index0, node_index0, h);
	fillGlobalStiffness(_K, Fq0dq1, node_index0, node_index1, h);

	if (!node0->onBorder)
	{
		MatrixXd fri_u0dq0 = 0.5*node0->normal.transpose()*Fq0dq0;
		MatrixXd fri_u0dq1 = 0.5*node0->normal.transpose()*Fq0dq1;
		
		fillGlobalStiffness(node0->u_friction, fri_u0dq0, cross_index0, node_index0, h);
		fillGlobalStiffness(node0->u_friction, fri_u0dq1, cross_index0, node_index1, h);
	}

	if (!node1->onBorder)
	{
		MatrixXd fri_u1dq0 = 0.5*node1->normal.transpose() * Fq1dq0;
		MatrixXd fri_u1dq1 = 0.5*node1->normal.transpose() *Fq1dq1;

		fillGlobalStiffness(node1->u_friction, fri_u1dq0, cross_index1, node_index0, h);
		fillGlobalStiffness(node1->u_friction, fri_u1dq1, cross_index1, node_index1, h);
	}

	/* Euler Part */
	if (!node0->onBorder)
	{
		double Fu0du0 = -Ks * w.norm()*w.norm() / delta_u;
		_K.push_back(T(cross_index0, cross_index0, -h*h*Fu0du0));
		
		MatrixXd Fq0du0 = Ks * w.norm() / delta_u * d;
		MatrixXd Fq1du0 = -Fq0du0;
		fillGlobalStiffness(_K, Fq0du0, node_index0, cross_index0, h);
		fillGlobalStiffness(_K, Fq1du0, node_index1, cross_index0, h);

		if (!node0->onBorder)
		{
			MatrixXd fri_u0du0 = 0.5*node0->normal.transpose()*Fq0du0;
			fillGlobalStiffness(node0->u_friction, fri_u0du0, cross_index0, cross_index0, h);
		}

		if (!node1->onBorder)
		{
			MatrixXd fri_u1du0 = 0.5*node1->normal.transpose()*Fq1du0;
			fillGlobalStiffness(node1->u_friction, fri_u1du0, cross_index1, cross_index0, h);
		}

		MatrixXd Fu0dq0 = Ks / delta_u * w.transpose();
		MatrixXd Fu0dq1 = -Fu0dq0;
		fillGlobalStiffness(_K, Fu0dq0, cross_index0, node_index0, h);
		fillGlobalStiffness(_K, Fu0dq1, cross_index0, node_index1, h);
	}

	if (!node1->onBorder)
	{		
		double Fu1du1 = -Ks * w.norm()*w.norm() / delta_u;
		_K.push_back(T(cross_index1, cross_index1, -h * h*Fu1du1));
		
		MatrixXd Fq1du1 = Ks * w.norm() / delta_u * d;
		MatrixXd Fq0du1 = -Fq1du1;
		fillGlobalStiffness(_K, Fq1du1, node_index1, cross_index1, h);
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

		MatrixXd Fu1dq1 = Ks / delta_u * w.transpose();
		MatrixXd Fu1dq0 = -Fu1dq1;
		fillGlobalStiffness(_K, Fu1dq1, cross_index1, node_index1, h);
		fillGlobalStiffness(_K, Fu1dq0, cross_index1, node_index0, h);
	}

	if (!node0->onBorder && !node1->onBorder)
	{
		double Fu1du0 = Ks * w.norm()*w.norm() / delta_u;
		double Fu0du1 = Fu1du0;

		_K.push_back(T(cross_index1, cross_index0, -h * h*Fu1du0));
		_K.push_back(T(cross_index0, cross_index1, -h * h*Fu0du1));
	}
}

void StretchSpring::solveV(vector<T>& _K, VectorXd& f, int nodes_size, double h)
{
	double l = (node1->position - node0->position).norm();
	double delta_v = abs(node1->v - node0->v);

	Vector3d w = (node1->position - node0->position) / delta_v;
	Vector3d d = (node1->position - node0->position); d.normalize();

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d * d.transpose();

	stretchEnergy = 0.5*Ks*delta_v*(w.norm() - 1)*(w.norm() - 1);

	int node_index0 = node0->node_index * 3;
	int node_index1 = node1->node_index * 3;

	int cross_index0 = nodes_size * 3 + node0->cross_index * 2;
	int cross_index1 = nodes_size * 3 + node1->cross_index * 2;

	//Compute and fill the force vector
	Vector3d Fq1 = -Ks * (w.norm() - 1)*d;
	Vector3d Fq0 = -Fq1;

	f.segment<3>(node_index0) += Fq0;
	f.segment<3>(node_index1) += Fq1;

	node0->compressForce -= 0.5*node0->normal.dot(Fq0);
	node1->compressForce -= 0.5*node1->normal.dot(Fq1);

	if (!node0->onBorder)
	{
		double Fv0 = -0.5*Ks*(w.norm()*w.norm() - 1);
		f(cross_index0 + 1) += Fv0;
	}

	if (!node1->onBorder)
	{
		double Fv1 = 0.5*Ks*(w.norm()*w.norm() - 1);
		f(cross_index1 + 1) += Fv1;
	}

	//Compute and fill the stiffness matrix
	/* Lagrange Part */
	MatrixXd Fq1dq1 = Ks / l * P - Ks / delta_v * I;
	MatrixXd Fq0dq0 = Fq1dq1;
	MatrixXd Fq1dq0 = -Fq1dq1;
	MatrixXd Fq0dq1 = Fq1dq0;

	fillGlobalStiffness(_K, Fq1dq1, node_index1, node_index1, h);
	fillGlobalStiffness(_K, Fq1dq0, node_index1, node_index0, h);
	fillGlobalStiffness(_K, Fq0dq0, node_index0, node_index0, h);
	fillGlobalStiffness(_K, Fq0dq1, node_index0, node_index1, h);

	if (!node0->onBorder)
	{
		MatrixXd fri_v0dq0 = -0.5*node0->normal.transpose()*Fq0dq0;
		MatrixXd fri_v0dq1 = -0.5*node0->normal.transpose()*Fq0dq1;

		fillGlobalStiffness(node0->v_friction, fri_v0dq0, cross_index0 + 1, node_index0, h);
		fillGlobalStiffness(node0->v_friction, fri_v0dq1, cross_index0 + 1, node_index1, h);
	}

	if (!node1->onBorder)
	{
		MatrixXd fri_v1dq0 = -0.5*node1->normal.transpose()*Fq1dq0;
		MatrixXd fri_v1dq1 = -0.5*node1->normal.transpose()*Fq1dq1;

		fillGlobalStiffness(node1->v_friction, fri_v1dq0, cross_index1 + 1, node_index0, h);
		fillGlobalStiffness(node1->v_friction, fri_v1dq1, cross_index1 + 1, node_index1, h);
	}

	/* Euler Part */
	if (!node0->onBorder)
	{
		double Fv0dv0 = -Ks * w.norm()*w.norm() / delta_v;
		_K.push_back(T(cross_index0 + 1, cross_index0 + 1, -h * h*Fv0dv0));

		MatrixXd Fq0dv0 = Ks * w.norm() / delta_v * d;
		MatrixXd Fq1dv0 = -Fq0dv0;
		fillGlobalStiffness(_K, Fq0dv0, node_index0, cross_index0 + 1, h);
		fillGlobalStiffness(_K, Fq1dv0, node_index1, cross_index0 + 1, h);

		if (!node0->onBorder)
		{
			MatrixXd fri_v0dv0 = -0.5*node0->normal.transpose()*Fq0dv0;
			fillGlobalStiffness(node0->v_friction, fri_v0dv0, cross_index0 + 1, cross_index0 + 1, h);
		}

		if (!node1->onBorder)
		{
			MatrixXd fri_v1dv0 = -0.5*node1->normal.transpose()*Fq1dv0;
			fillGlobalStiffness(node1->v_friction, fri_v1dv0, cross_index1 + 1, cross_index0 + 1, h);
		}

		MatrixXd Fv0dq0 = Ks / delta_v * w.transpose();
		MatrixXd Fv0dq1 = -Fv0dq0;
		fillGlobalStiffness(_K, Fv0dq0, cross_index0 + 1, node_index0, h);
		fillGlobalStiffness(_K, Fv0dq1, cross_index0 + 1, node_index1, h);
	}

	if (!node1->onBorder)
	{
		double Fv1dv1 = -Ks * w.norm()*w.norm() / delta_v;
		_K.push_back(T(cross_index1 + 1, cross_index1 + 1, -h * h*Fv1dv1));

		MatrixXd Fq1dv1 = Ks * w.norm() / delta_v * d;
		MatrixXd Fq0dv1 = -Fq1dv1;
		fillGlobalStiffness(_K, Fq1dv1, node_index1, cross_index1 + 1, h);
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

		MatrixXd Fv1dq1 = Ks / delta_v * w.transpose();
		MatrixXd Fv1dq0 = -Fv1dq1;
		fillGlobalStiffness(_K, Fv1dq1, cross_index1 + 1, node_index1, h);
		fillGlobalStiffness(_K, Fv1dq0, cross_index1 + 1, node_index0, h);
	}

	if (!node0->onBorder && !node1->onBorder)
	{
		double Fv1dv0 = Ks * w.norm()*w.norm() / delta_v;;
		double Fv0dv1 = Fv1dv0;

		_K.push_back(T(cross_index1, cross_index0, -h * h*Fv1dv0));
		_K.push_back(T(cross_index0, cross_index1, -h * h*Fv0dv1));
	}
}

void StretchSpring::draw()
{
	if (springType == Weft)
	{
		glColor3d(0.9, 0.9, 0.9);
		glLineWidth(2.0f);
		glBegin(GL_LINES);
		Vector3d pos0, pos1;
		if (node0->neighbor[NodeLocation::Up] && node0->neighbor[NodeLocation::Down] && node0->neighbor[NodeLocation::Right] && node0->neighbor[NodeLocation::Left])
		{
			pos0 = node0->position + node0->normal*R/2.0;
		}
		else
		{
			pos0 = node0->position;
		}

		if (node1->neighbor[NodeLocation::Up] && node1->neighbor[NodeLocation::Down] && node1->neighbor[NodeLocation::Right] && node1->neighbor[NodeLocation::Left])
		{
			pos1 = node1->position + node1->normal*R/2.0;
		}
		else
		{
			pos1 = node1->position;
		}
		glVertex3d(pos0.x()*1e3, pos0.y()*1e3, pos0.z()*1e3);
		glVertex3d(pos1.x()*1e3, pos1.y()*1e3, pos1.z()*1e3);
		glEnd();
	}
	else if (springType == Warp)
	{
		glColor3d(0.9, 0.6, 0.3);
		glLineWidth(2.0f);
		glBegin(GL_LINES);
		Vector3d pos0, pos1;
		if (node0->neighbor[NodeLocation::Up] && node0->neighbor[NodeLocation::Down] && node0->neighbor[NodeLocation::Right] && node0->neighbor[NodeLocation::Left])
		{
			pos0 = node0->position - node0->normal*R/2.0;
		}
		else
		{
			pos0 = node0->position;
		}

		if (node1->neighbor[NodeLocation::Up] && node1->neighbor[NodeLocation::Down] && node1->neighbor[NodeLocation::Right] && node1->neighbor[NodeLocation::Left])
		{
			pos1 = node1->position - node1->normal*R/2.0;
		}
		else
		{
			pos1 = node1->position;
		}		
		glVertex3d(pos0.x()*1e3, pos0.y()*1e3, pos0.z()*1e3);
		glVertex3d(pos1.x()*1e3, pos1.y()*1e3, pos1.z()*1e3);
		glEnd();
	}
}
