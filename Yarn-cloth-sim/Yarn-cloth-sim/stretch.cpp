#include "stretch.h"
#include <freeglut.h>

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

void StretchSpring::solve(vector<T>& _K, VectorXd& f)
{
	if (springType == Weft) solveV(_K, f);
	
	else if (springType == Warp) solveU(_K, f);

	return;
}

void StretchSpring::solveU(vector<T>& _K, VectorXd& f)
{
	double l = (node1->position - node0->position).norm();
	double delta_u = abs(node1->u - node0->u);

	Vector3d w = (node1->position - node0->position) / delta_u;
	Vector3d d = node1->position - node0->position; d.normalize();
	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d * d.transpose();

	double V = 0.5*Ks*delta_u*(w.norm() - 1)*(w.norm() - 1);

	int index0 = node0->index * 5;
	int index1 = node1->index * 5;

	/* Compute and fill the force vector */
	Vector3d Fq1 = -Ks * (w.norm() - 1)*d;
	Vector3d Fq0 = -Fq1;

	double Fu1 = 0.5*Ks*(w.norm()*w.norm() - 1);
	double Fu0 = -Fu1;
	
	f.segment<3>(index0) += Fq0;
	f(index0 + 3) += Fu0;

	f.segment<3>(index1) += Fq1;
	f(index1 + 3) += Fu1;

	/* Compute and fill the stiffness matrix */
	//The local stiffness matrix should be 10*10
	Matrix3d Fq1dq1 = Ks / l * P - Ks / delta_u * I;
	Matrix3d Fq0dq0 = Fq1dq1;
	Matrix3d Fq1dq0 = -Fq1dq1;
	Matrix3d Fq0dq1 = Fq1dq0;

	double Fu1du1 = -Ks * w.norm()*w.norm() / delta_u;
	double Fu0du0 = Fu1du1;
	double Fu1du0 = -Fu1du1;
	double Fu0du1 = Fu1du0;

	Vector3d Fq1du1 = Ks * w.norm() / delta_u * d;
	Vector3d Fq0du0 = Fq1du1;
	Vector3d Fq1du0 = -Fq1du1;
	Vector3d Fq0du1 = Fq1du0;

	Vector3d Fu1dq1 = Ks / delta_u * w.transpose();
	Vector3d Fu0dq0 = Fu1dq1;
	Vector3d Fu1dq0 = -Fu1dq1;
	Vector3d Fu0dq1 = Fu1dq0;

	MatrixXd f0q0;
	f0q0.resize(5, 5);
	f0q0.setZero();
	f0q0.block<3, 3>(0, 0) = Fq0dq0;
	f0q0.block<3, 1>(0, 3) = Fq0du0;
	f0q0.block<1, 3>(3, 0) = Fu0dq0;
	f0q0(3, 3) = Fu0du0;
	fillBlock(_K, f0q0, index0, index0);

	MatrixXd f0q1;
	f0q1.resize(5, 5);
	f0q1.setZero();
	f0q1.block<3, 3>(0, 0) = Fq0dq1;
	f0q1.block<3, 1>(0, 3) = Fq0du1;
	f0q1.block<1, 3>(3, 0) = Fu0dq1;
	f0q1(3, 3) = Fu0du1;
	fillBlock(_K, f0q1, index0, index1);

	MatrixXd f1q0;
	f1q0.resize(5, 5);
	f1q0.setZero();
	f1q0.block<3, 3>(0, 0) = Fq1dq0;
	f1q0.block<3, 1>(0, 3) = Fq1du0;
	f1q0.block<1, 3>(3, 0) = Fu1dq0;
	f1q0(3, 3) = Fu1du0;
	fillBlock(_K, f1q0, index1, index0);

	MatrixXd f1q1;
	f1q1.resize(5, 5);
	f1q1.setZero();
	f1q1.block<3, 3>(0, 0) = Fq1dq1;
	f1q1.block<3, 1>(0, 3) = Fq1du1;
	f1q1.block<1, 3>(3, 0) = Fu1dq1;
	f1q1(3, 3) = Fu1du1;
	fillBlock(_K, f1q1, index1, index1);
}

void StretchSpring::solveV(vector<T>& _K, VectorXd& f)
{
	double l = (node1->position - node0->position).norm();
	double delta_v = abs(node1->v - node0->v);
	
	Vector3d w = (node1->position - node0->position) / delta_v;
	Vector3d d = (node1->position - node0->position); d.normalize();

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d * d.transpose();

	double V = 0.5*Ks*delta_v*(w.norm() - 1)*(w.norm() - 1);

	int index0 = node0->index * 5;
	int index1 = node1->index * 5;

	//Compute and fill the force vector
	Vector3d Fq1 = -Ks * (w.norm() - 1)*d;
	Vector3d Fq0 = -Fq1;

	double Fv1 = 0.5*Ks*(w.norm()*w.norm() - 1);
	double Fv0 = -Fv1;

	f.segment<3>(index0) += Fq0;
	f(index0 + 4) += Fv0;

	f.segment<3>(index1) += Fq1;
	f(index1 + 4) += Fv1;

	//Compute and fill the stiffness matrix
	//The local stiffness matrix should be 10*10
	Matrix3d Fq1dq1 = Ks / l * P - Ks / delta_v * I;
	Matrix3d Fq0dq0 = Fq1dq1;
	Matrix3d Fq1dq0 = -Fq1dq1;
	Matrix3d Fq0dq1 = Fq1dq0;

	double Fv1dv1 = -Ks * w.norm()*w.norm() / delta_v;
	double Fv0dv0 = Fv1dv1;
	double Fv1dv0 = -Fv1dv1;
	double Fv0dv1 = Fv1dv0;

	Vector3d Fq1dv1 = Ks * w.norm() / delta_v * d;
	Vector3d Fq0dv0 = Fq1dv1;
	Vector3d Fq1dv0 = -Fq1dv1;
	Vector3d Fq0dv1 = Fq1dv0;

	Vector3d Fv1dq1 = Ks / delta_v * w.transpose();
	Vector3d Fv0dq0 = Fv1dq1;
	Vector3d Fv1dq0 = -Fv1dq1;
	Vector3d Fv0dq1 = Fv1dq0;

	MatrixXd f0q0;
	f0q0.resize(5, 5);
	f0q0.setZero();
	f0q0.block<3, 3>(0, 0) = Fq0dq0;
	f0q0.block<3, 1>(0, 4) = Fq0dv0;
	f0q0.block<1, 3>(4, 0) = Fv0dq0;
	f0q0(4, 4) = Fv0dv0;
	fillBlock(_K, f0q0, index0, index0);

	MatrixXd f0q1;
	f0q1.resize(5, 5);
	f0q1.setZero();
	f0q1.block<3, 3>(0, 0) = Fq0dq1;
	f0q1.block<3, 1>(0, 4) = Fq0dv1;
	f0q1.block<1, 3>(4, 0) = Fv0dq1;
	f0q1(4, 4) = Fv0dv1;
	fillBlock(_K, f0q1, index0, index1);

	MatrixXd f1q0;
	f1q0.resize(5, 5);
	f1q0.setZero();
	f1q0.block<3, 3>(0, 0) = Fq1dq0;
	f1q0.block<3, 1>(0, 4) = Fq1dv0;
	f1q0.block<1, 3>(4, 0) = Fv1dq0;
	f1q0(4, 4) = Fv1dv0;
	fillBlock(_K, f1q0, index1, index0);

	MatrixXd f1q1;
	f1q1.resize(5, 5);
	f1q1.setZero();
	f1q1.block<3, 3>(0, 0) = Fq1dq1;
	f1q1.block<3, 1>(0, 4) = Fq1dv1;
	f1q1.block<1, 3>(4, 0) = Fv1dq1;
	f1q1(4, 4) = Fv1dv1;
	fillBlock(_K, f1q1, index1, index1);
}

void StretchSpring::draw()
{
	if (springType == Weft)
	{
		if (node0->whichUp == Weft)
		{
			glColor3d(0.9, 0.9, 0.9);
			glLineWidth(2.0f);
			glBegin(GL_LINES);
			Vector3d pos0 = node0->position + node0->getNormal()*R;
			Vector3d pos1 = node1->position - node1->getNormal()*R;
			glVertex3d(pos0.x(), pos0.y(), pos0.z());
			glVertex3d(pos1.x(), pos1.y(), pos1.z());
			glEnd();
		}
		else if (node1->whichUp==Weft)
		{
			glColor3d(0.9, 0.9, 0.9);
			glLineWidth(2.0f);
			glBegin(GL_LINES);
			Vector3d pos0 = node0->position - node0->getNormal()*R;
			Vector3d pos1 = node1->position + node1->getNormal()*R;
			glVertex3d(pos0.x(), pos0.y(), pos0.z());
			glVertex3d(pos1.x(), pos1.y(), pos1.z());
			glEnd();
		}
	}
	else if (springType==Warp)
	{
		if (node0->whichUp == Warp)
		{
			glColor3d(0.9, 0.9, 0.9);
			glLineWidth(2.0f);
			glBegin(GL_LINES);
			Vector3d pos0 = node0->position + node0->getNormal()*R;
			Vector3d pos1 = node1->position - node1->getNormal()*R;
			glVertex3d(pos0.x(), pos0.y(), pos0.z());
			glVertex3d(pos1.x(), pos1.y(), pos1.z());
			glEnd();
		}
		else if (node1->whichUp == Warp)
		{
			glColor3d(0.9, 0.9, 0.9);
			glLineWidth(2.0f);
			glBegin(GL_LINES);
			Vector3d pos0 = node0->position - node0->getNormal()*R;
			Vector3d pos1 = node1->position + node1->getNormal()*R;
			glVertex3d(pos0.x(), pos0.y(), pos0.z());
			glVertex3d(pos1.x(), pos1.y(), pos1.z());
			glEnd();
		}
	}
}
