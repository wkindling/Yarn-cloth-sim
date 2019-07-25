#pragma once
#ifndef NODE_H
#define NODE_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "utility.h"

//Left & Right for Warp, Up & Down for Weft
enum NodeLocation
{
	Left,
	Up,
	Right,
	Down
};//ClockWise

enum YarnType // Warp for horizontal, Weft for vertical
{
	//Normal,
	Warp,
	Weft,
};

//Currently we ignore the normal nodes in the cloth
class Node
{
public:
	Node(Eigen::Vector3d pos, double _u, double _v);
	virtual ~Node();

	Eigen::Vector3d getNormal();

public:
	//Normal Node
	Eigen::Vector3d position;
	Eigen::Vector3d velocity;
	Eigen::Vector3d normal;

	Node* neighbor[4];

	//Cross Node
	double u, v; // u for warp, v for weft. 
	double anchorU, anchorV; // Initial u,v
	
	Eigen::Vector2d velocityUV;

	YarnType whichUp;

	int index; // Record global matrix

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif