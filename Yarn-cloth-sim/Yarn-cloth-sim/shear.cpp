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

void ShearSpring::solve(vector<T>& _K, VectorXd& f)
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

	shearEnergy = 0.5*Kx*L*(phi - M_PI / 2.0)*(phi - M_PI / 2.0);

	int index0 = node0->index * 5;
	int index1 = node1->index * 5;
	int index3 = node3->index * 5;

	//Compute and fill the force vector
	//Forces on warp and weft coordinates are all zero.
	Vector3d Fq1 = Kx * L*(phi - M_PI / 2.0) / l1 / sin(phi)*P1*d3;
	Vector3d Fq3 = Kx * L*(phi - M_PI / 2.0) / l3 / sin(phi)*P3*d1;
	Vector3d Fq0 = -(Fq1 + Fq3);

	f.segment<3>(index0) += Fq0;
	f.segment<3>(index1) += Fq1;
	f.segment<3>(index3) += Fq3;
	
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

	MatrixXd f0q0;
	f0q0.resize(5, 5);
	f0q0.setZero();
	f0q0.block<3, 3>(0, 0) = Fq0dq0;
	fillBlock(_K, f0q0, index0, index0);

	MatrixXd f0q1;
	f0q1.resize(5, 5);
	f0q1.setZero();
	f0q1.block<3, 3>(0, 0) = Fq0dq1;
	fillBlock(_K, f0q1, index0, index1);

	MatrixXd f0q3;
	f0q3.resize(5, 5);
	f0q3.setZero();
	f0q3.block<3, 3>(0, 0) = Fq0dq3;
	fillBlock(_K, f0q3, index0, index3);

	MatrixXd f1q0;
	f1q0.resize(5, 5);
	f1q0.setZero();
	f1q0.block<3, 3>(0, 0) = Fq1dq0;
	fillBlock(_K, f1q0, index1, index0);

	MatrixXd f1q1;
	f1q1.resize(5, 5);
	f1q1.setZero();
	f1q1.block<3, 3>(0, 0) = Fq1dq1;
	fillBlock(_K, f1q1, index1, index1);

	MatrixXd f1q3;
	f1q3.resize(5, 5);
	f1q3.setZero();
	f1q3.block<3, 3>(0, 0) = Fq1dq3;
	fillBlock(_K, f1q3, index1, index3);

	MatrixXd f3q0;
	f3q0.resize(5, 5);
	f3q0.setZero();
	f3q0.block<3, 3>(0, 0) = Fq3dq0;
	fillBlock(_K, f3q0, index3, index0);

	MatrixXd f3q1;
	f3q1.resize(5, 5);
	f3q1.setZero();
	f3q1.block<3, 3>(0, 0) = Fq3dq1;
	fillBlock(_K, f3q1, index3, index1);

	MatrixXd f3q3;
	f3q3.resize(5, 5);
	f3q3.setZero();
	f3q3.block<3, 3>(0, 0) = Fq3dq3;
	fillBlock(_K, f3q3, index3, index3);
}