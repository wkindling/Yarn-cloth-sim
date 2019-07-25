#include "node.h"
#include <iostream>

using namespace Eigen;
using namespace std;

Node::Node(Vector3d pos, double _u, double _v)
{
	position = pos;
	velocity.setZero();

	u = _u;
	v = _v;
	anchorU = _u;
	anchorV = _v;
	velocityUV.setZero();

	normal.setZero();
	compressForce = 0;
	//whichUp and index assigned externally
}

Node::~Node() {}


// Use SVD to get the vertex normal according to neighbors
// Ensure that normal points from warp to weft
Vector3d Node::getNormal()
{
	/* SVD */
	int neighbor_count = 0;
	Vector3d mean=Vector3d::Zero();
	mean += this->position;
	for (int i = 0; i < 4; i++)
	{
		if (neighbor[i] != NULL)
		{
			mean += neighbor[i]->position;
			neighbor_count++;
		}
	}
	mean = mean / (double)(1 + neighbor_count);
	MatrixXd A;
	A.resize(neighbor_count + 1, 3);

	A.row(0) = (this->position - mean).transpose();
	int row_index = 1;
	for (int i = 0; i < 4; i++)
	{
		if (neighbor[i] != NULL)
		{
			A.row(row_index) = (neighbor[i]->position - mean).transpose();
			row_index++;
		}
	}
	
	JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);
	MatrixXd V = svd.matrixV();
	Vector3d nor = V.col(2);
	nor.normalize();
	
	/* Ensure that normal points from warp to weft */
	Vector3d cloth_up;
	if (neighbor[NodeLocation::Up] && neighbor[NodeLocation::Right])
	{
		Vector3d up = neighbor[NodeLocation::Up]->position - this->position;
		Vector3d right = neighbor[NodeLocation::Right]->position - this->position;
		cloth_up = right.cross(up);
		cloth_up.normalize();
	}
	else if (neighbor[NodeLocation::Right] && neighbor[NodeLocation::Down])
	{
		Vector3d right = neighbor[NodeLocation::Right]->position - this->position;
		Vector3d down = neighbor[NodeLocation::Down]->position - this->position;
		cloth_up = down.cross(right);
		cloth_up.normalize();
	}
	else if (neighbor[NodeLocation::Down] && neighbor[NodeLocation::Left])
	{
		Vector3d down = neighbor[NodeLocation::Down]->position - this->position;
		Vector3d left = neighbor[NodeLocation::Left]->position - this->position;
		cloth_up = left.cross(down);
		cloth_up.normalize();
	}
	else if (neighbor[NodeLocation::Left] && neighbor[NodeLocation::Up])
	{
		Vector3d up = neighbor[NodeLocation::Up]->position - this->position;
		Vector3d left = neighbor[NodeLocation::Left]->position - this->position;
		cloth_up = up.cross(left);
		cloth_up.normalize();
	}

	if (whichUp == Weft) 
	{
		nor = nor.dot(cloth_up) > 0 ? nor : -nor;
	}
	else if (whichUp == Warp)
	{
		nor = nor.dot(cloth_up) < 0 ? nor : -nor;
	}

	normal = nor;

	return normal;
}

/* Compute friction force */
// Ignore damping force yet ... Because damping force may have some influence on global motion equation
// Currently I ignore the computation of Jacobians for friction because it is too complicated
double Node::getFriction(double mu, double Kf, vector<T>& _K, VectorXd& f)
{
	double friction_u = 0, friction_v = 0;
	double limit = mu * compressForce;

	//Maybe I need to fill the stiffness matrix again due to friction
	/* Get friction in warp direction */
	double hat_fu = -Kf * (u - anchorU);
	
	if (hat_fu <= limit)
	{
		friction_u = hat_fu;
	}
	else
	{
		friction_u = -sign(u - anchorU)*mu*compressForce;
	}

	f(index * 5 + 3) += friction_u;

	/* Get friction in weft direction */
	double hat_fv = -Kf * (v - anchorV);
	
	if (hat_fv <= limit)
	{
		friction_v = hat_fv;
	}
	else
	{
		friction_v = -sign(v - anchorV)*mu*compressForce;
	}

	f(index * 5 + 4) += friction_v;
}