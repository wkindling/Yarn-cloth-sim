#include "shear.h"

using namespace std;
using namespace Eigen;

//n0-n1 for warp 'u', and n0-n3 for weft 'v'
ShearSpring::ShearSpring(Node* n0, Node* n1, Node* n3, double _S, double _R, double _L)
{
	node0 = n0;
	node1 = n1;
	node3 = n3;
	S = _S;
	R = _R;
	L = _L;
	
	Kx = S * R*R;
	shearEnergy = 0;
}

ShearSpring::~ShearSpring() {}

void ShearSpring::solve()
{
	double l1 = (node1->position - node0->position).norm();
	double l3 = (node3->position - node0->position).norm();

	Vector3d d1 = node1->position - node0->position; d1.normalize();
	Vector3d d3 = node3->position - node0->position; d3.normalize();

	double delta_u = abs(node1->u - node0->u);
	double delta_v = abs(node3->v - node0->v);

	Vector3d w1 = (node1->position - node0->position) / delta_u;
	Vector3d w3 = (node3->position - node0->position) / delta_v;

	Matrix3d I = Matrix3d::Identity();
	Matrix3d P1 = I - d1 * d1.transpose();
	Matrix3d P3 = I - d3 * d3.transpose();

	double phi = acos(d1.dot(d3));

	double V = 0.5*Kx*L*(phi - M_PI / 2.0)*(phi - M_PI / 2.0);

	//Compute and fill the force vector
	//Forces on warp and weft coordinates are all zero.
	Vector3d Fq1 = Kx * L*(phi - M_PI / 2.0) / l1 / sin(phi)*P1*d3;
	Vector3d Fq3 = Kx * L*(phi - M_PI / 2.0) / l3 / sin(phi)*P3*d1;
	Vector3d Fq0 = -(Fq1 + Fq3);

	//Compute and fill the stiffness matrix
	//The local stiffness matrix should be 15*15, but the derivatives of u and v are zero 
	Matrix3d Fq1dq1 = Kx * L / l1 / l1 / sin(phi)*((phi - M_PI / 2.0)*(-P1 * d3*d1.transpose() + cos(phi) / sin(phi) / sin(phi)*P1*d3*d3.transpose()*P1 - cos(phi)*P1 - d1 * d3.transpose()*P1) - 1 / sin(phi)*P1*d3*d3.transpose()*P1);
	Matrix3d Fq1dq3 = Kx * L / l3 / l1 / sin(phi)*((phi - M_PI / 2.0)*(cos(phi) / sin(phi) / sin(phi)*P1*d3*d1.transpose() + P1) - 1 / sin(phi)*P1*d3*d1.transpose())*P3;
	Matrix3d Fq3dq1 = Kx * L / l1 / l3 / sin(phi)*((phi - M_PI / 2.0)*(cos(phi) / sin(phi) / sin(phi)*P3*d1*d3.transpose() + P3) - 1 / sin(phi)*P3*d1*d3.transpose())*P1;
	Matrix3d Fq3dq3 = Kx * L / l3 / l3 / sin(phi)*((phi - M_PI / 2.0)*(-P3 * d1*d3.transpose() + cos(phi) / sin(phi) / sin(phi)*P3*d1*d1.transpose()*P3 - cos(phi)*P3 - d3 * d1.transpose()*P3) - 1 / sin(phi)*P3*d1*d1.transpose()*P3);

	Matrix3d Fq1dq0 = -(Fq1dq1 + Fq1dq3);
	Matrix3d Fq3dq0 = -(Fq3dq1 + Fq3dq3);
	Matrix3d Fq0dq1 = -(Fq1dq1 + Fq3dq1);
	Matrix3d Fq0dq3 = -(Fq1dq3 + Fq3dq3);
	Matrix3d Fq0dq0 = -(Fq1dq0 + Fq3dq0);

	/* TODO: Fill the block just like EoL*/




}