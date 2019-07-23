#include "stretch.h"
#include <freeglut.h>

using namespace std;
using namespace Eigen;

StretchSpring::StretchSpring(Node* n0, Node* n1, double _Y, double _R, YarnType type)
{
	node0 = n0;
	node1 = n1;
	Y = _Y;
	R = _R;

	Ks = Y * M_PI*R*R;
	stretchEnergy = 0;
	springType = type;
}

StretchSpring::~StretchSpring() {}

void StretchSpring::solve()
{
	if (springType == Weft) solveV();
	
	else if (springType == Warp) solveU();

	return;
}

void StretchSpring::solveU()
{
	double l = (node1->position - node0->position).norm();
	double delta_u = abs(node1->u - node0->u);

	Vector3d w = (node1->position - node0->position) / delta_u;
	Vector3d d = node1->position - node0->position; d.normalize();
	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d * d.transpose();

	//Compute and fill the force vector
	Vector3d Fq1 = -Ks * (w.norm() - 1)*d;
	Vector3d Fq0 = -Fq1;

	double Fu1 = 0.5*Ks*(w.norm()*w.norm() - 1);
	double Fu0 = -Fu1;

	//Compute and fill the stiffness matrix
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

	/* TODO: Fill the block just like EoL*/
	
}

void StretchSpring::solveV()
{
	double l = (node1->position - node0->position).norm();
	double delta_v = abs(node1->v - node0->v);
	
	Vector3d w = (node1->position - node0->position) / delta_v;
	Vector3d d = (node1->position - node0->position); d.normalize();

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P = I - d * d.transpose();

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
