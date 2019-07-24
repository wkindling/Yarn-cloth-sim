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

	f(index0) = Fq0.x();
	f(index0 + 1) = Fq0.y();
	f(index0 + 2) = Fq0.z();
	f(index0 + 3) = Fu0;
	f(index0 + 4) = 0;
	
	f(index1) = Fq1.x();
	f(index1 + 1) = Fq1.y();
	f(index1 + 2) = Fq1.z();
	f(index1 + 3) = Fu1;
	f(index1 + 4) = 0;

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

	/* Fill the block just like EoL*/
	//Fill Fq0 related blocks
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_K.push_back(T(index0 + i, index0 + j, Fq0dq0(i, j))); //Fq0 dq0
			_K.push_back(T(index0 + i, index1 + j, Fq0dq1(i, j))); //Fq0 dq1
		}

		_K.push_back(T(index0 + i, index0 + 3, Fq0du0(i))); //Fq0 du0
		_K.push_back(T(index0 + i, index0 + 4, 0)); //Fq0 dv0
		_K.push_back(T(index0 + i, index1 + 3, Fq0du1(i))); //Fq0 du1
		_K.push_back(T(index0 + i, index1 + 4, 0)); //Fq0 dv1
	}

	//Fill Fq1 related blocks
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			_K.push_back(T(index1 + i, index0 + j, Fq1dq0(i, j))); //Fq1 q0
			_K.push_back(T(index1 + i, index1 + j, Fq1dq1(i, j))); //Fq1 dq1
		}

		_K.push_back(T(index1 + i, index0 + 3, Fq1du0(i))); //Fq1 du0
		_K.push_back(T(index1 + i, index0 + 4, 0)); //Fq1 dv0
		_K.push_back(T(index1 + i, index1 + 3, Fq1du1(i))); //Fq1 du1
		_K.push_back(T(index1 + i, index1 + 4, 0)); //Fq1 dv1
	}

	for (int i = 0; i < 3; i++)
	{
		_K.push_back(T(index0 + 3, index0 + i, Fu0dq0(i))); //Fu0 dq0
		_K.push_back(T(index0 + 4, index0 + i, 0)); //Fv0 dq0
		_K.push_back(T(index0 + 3, index1 + i, Fu0dq1(i))); //Fu0 dq1
		_K.push_back(T(index0 + 4, index1 + i, 0)); //Fv0 dq1
	}

	_K.push_back(T(index0 + 3, index0 + 3, Fu0du0)); //Fu0 du0
	_K.push_back(T(index0 + 3, index0 + 4, 0)); //Fu0 dv0
	_K.push_back(T(index0 + 4, index0 + 3, 0)); //Fv0 du0
	_K.push_back(T(index0 + 4, index0 + 4, 0)); //Fv0 dv0

	_K.push_back(T(index0 + 3, index1 + 3, Fu0du1));//Fu0 du1
	_K.push_back(T(index0 + 3, index1 + 4, 0));//Fu0 dv1
	_K.push_back(T(index0 + 4, index1 + 3, 0)); //Fv0 du1
	_K.push_back(T(index0 + 4, index1 + 4, 0)); //Fv0 dv1

	for (int i = 0; i < 3; i++)
	{
		_K.push_back(T(index1 + 3, index0 + i, Fu1dq0(i))); //Fu1 dq0
		_K.push_back(T(index1 + 4, index0 + i, 0)); //Fv1 dq0
		_K.push_back(T(index1 + 3, index1 + i, Fu1dq1(i))); //Fu1 dq1
		_K.push_back(T(index1 + 4, index1 + i, 0)); //Fv1 dq1
	}

	_K.push_back(T(index1 + 3, index0 + 3, Fu1du0)); //Fu1 du0
	_K.push_back(T(index1 + 3, index0 + 4, 0)); //Fu1 dv0
	_K.push_back(T(index1 + 4, index0 + 3, 0)); //Fv1 du0
	_K.push_back(T(index1 + 4, index0 + 4, 0)); //Fv1 dv0

	_K.push_back(T(index1 + 3, index1 + 3, Fu1du1)); // Fu1 du1
	_K.push_back(T(index1 + 3, index1 + 4, 0)); // Fu1 dv1 
	_K.push_back(T(index1 + 4, index1 + 3, 0)); // Fv1 du1
	_K.push_back(T(index1 + 4, index1 + 4, 0)); //Fv1 dv1

	
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

	//Compute and fill the force vector
	Vector3d Fq1 = -Ks * (w.norm() - 1)*d;
	Vector3d Fq0 = -Fq1;

	double Fv1 = 0.5*Ks*(w.norm()*w.norm() - 1);
	double Fv0 = -Fv1;



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

	/*TODO: Fill the block just like EoL*/








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
			glVertex3d(node0->position.x(), node0->position.y(), node0->position.z() + R);
			glVertex3d(node1->position.x(), node1->position.y(), node1->position.z() - R);
			glEnd();
		}
		else if (node1->whichUp==Weft)
		{
			glColor3d(0.9, 0.9, 0.9);
			glLineWidth(2.0f);
			glBegin(GL_LINES);
			glVertex3d(node0->position.x(), node0->position.y(), node0->position.z() - R);
			glVertex3d(node1->position.x(), node1->position.y(), node1->position.z() + R);
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
			glVertex3d(node0->position.x(), node0->position.y(), node0->position.z() + R);
			glVertex3d(node1->position.x(), node1->position.y(), node1->position.z() - R);
			glEnd();
		}
		else if (node1->whichUp == Warp)
		{
			glColor3d(0.9, 0.9, 0.9);
			glLineWidth(2.0f);
			glBegin(GL_LINES);
			glVertex3d(node0->position.x(), node0->position.y(), node0->position.z() - R);
			glVertex3d(node1->position.x(), node1->position.y(), node1->position.z() + R);
			glEnd();
		}
	}
}
