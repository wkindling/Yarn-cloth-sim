#pragma once
#ifndef NODE_H
#define NODE_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#define M_PI 3.1415926535

//Left & Right for Warp, Up & Down for Weft
enum NodeLocation
{
	Left,
	Up,
	Right,
	Down
};//ClockWise

enum YarnType // Whcih yarn is up, normal for single yarn
{
	Normal,
	Warp,
	Weft,
};

class Node
{
public:
	Node(Eigen::Vector3d pos, double _u, double _v);
	virtual ~Node();

public:
	//Normal Node
	Eigen::Vector3d position;
	Eigen::Vector3d velocity;

	Node* neighbor[4];

	//Cross Node
	double u, v; // u for warp, v for weft. 
	double anchorU, anchorV; // Initial u,v
	
	Eigen::Vector2d velocityUV;

	YarnType whichUp;

	int cross_index;
	int normal_index; // Record global matrix

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


#endif