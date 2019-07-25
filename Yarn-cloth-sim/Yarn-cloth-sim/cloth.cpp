#include "cloth.h"
#include "node.h"
#include "stretch.h"
#include "bend.h"
#include "parallel_contact.h"
#include "shear.h"

#include <iostream>

using namespace std;
using namespace Eigen;

Cloth::Cloth(int w, int h, double _R, double _L)
{
	width = w;
	height = h;
	R = _R;
	L = _L;

	int cindex = 0;
	int nindex = 0;

	cout << "Cloth Width : " << width << endl;
	cout << "Cloth Height : " << height << endl;
	cout << "Yarn Radius : " << R << endl;
	cout << "Inter-Yarn Distance : " << L << endl;
	
	build();

}

Cloth::~Cloth() {}

void Cloth::build()
{
	//Create node list
	for (int j = 0; j < height; j++) //u for warp (x)
	{
		for (int i = 0; i < width; i++) // v for weft (y)
		{
			Node* node = new Node(Vector3d(i*L, j*L, 0), i*L, j*L);

			if ((i + j) % 2 == 0)
			{
				node->whichUp = Warp;
			}
			else node->whichUp = Weft;

			node->index = nodes.size(); //Record index for later global matrix
			nodes.push_back(node);
		}
	}

	for (int j = 0; j < height; j++) //Record the neighbors of each node
	{
		for (int i = 0; i < width; i++)
		{
			Node* node = nodes[j*width + i];
			node->neighbor[NodeLocation::Left] = i >= 1 ? nodes[j*width + i - 1] : NULL;
			node->neighbor[NodeLocation::Right] = i < width - 1 ? nodes[j*width + i + 1] : NULL;
			node->neighbor[NodeLocation::Up] = j >= 1 ? nodes[(j - 1)*width + i] : NULL;
			node->neighbor[NodeLocation::Down] = j < height - 1 ? nodes[(j + 1)*width + i] : NULL;
		}
	}

	for (int i = 0; i < nodes.size(); i++)
	{
		nodes[i]->getNormal();
	}

	/*Add 4 kinds of Springs */
	//Add stretch spring, avoid duplicate
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			Node* p = nodes[j*width + i];
			Node* right = p->neighbor[NodeLocation::Right];
			Node* forward = p->neighbor[NodeLocation::Down];
			if (right)
			{
				StretchSpring* warpStretch = new StretchSpring(p, right, Y, R, Warp);
				stretch_springs.push_back(warpStretch);
			}
			if (forward)
			{
				StretchSpring* weftStretch = new StretchSpring(p, forward, Y, R, Weft);
				stretch_springs.push_back(weftStretch);
			}
		}
	}

	//Add Bending Springs, no need to avoid duplicate
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			Node* p = nodes[j*width + i];
			if (p->neighbor[NodeLocation::Left] && p->neighbor[NodeLocation::Right])
			{
				BendSpring* warpBend = new BendSpring(p, p->neighbor[NodeLocation::Right], p->neighbor[NodeLocation::Left], B, R, Warp);
				bend_springs.push_back(warpBend);
			}
			if (p->neighbor[NodeLocation::Up] && p->neighbor[NodeLocation::Down])
			{
				BendSpring* weftBend = new BendSpring(p, p->neighbor[NodeLocation::Up], p->neighbor[NodeLocation::Down], B, R, Weft);
				bend_springs.push_back(weftBend);
			}
		}
	}

	//Add Shear Springs, no need to avoid duplicate
	//We need to ensure that n0-n3 is the weft yarn, and n0-n1 is the wap yarn
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			Node* p = nodes[j*width + i];
			if (p->neighbor[NodeLocation::Left] && p->neighbor[NodeLocation::Up])
			{
				ShearSpring* leftUpShear = new ShearSpring(p, p->neighbor[NodeLocation::Left], p->neighbor[NodeLocation::Up], S, R, L);
				shear_springs.push_back(leftUpShear);
			}

			if (p->neighbor[NodeLocation::Up] && p->neighbor[NodeLocation::Right])
			{
				ShearSpring* rightUpShear = new ShearSpring(p, p->neighbor[NodeLocation::Right], p->neighbor[NodeLocation::Up], S, R, L);
				shear_springs.push_back(rightUpShear);
			}

			if (p->neighbor[NodeLocation::Right] && p->neighbor[NodeLocation::Down])
			{
				ShearSpring* rightDownShear = new ShearSpring(p, p->neighbor[NodeLocation::Right], p->neighbor[NodeLocation::Down], S, R, L);
				shear_springs.push_back(rightDownShear);
			}

			if (p->neighbor[NodeLocation::Down] && p->neighbor[NodeLocation::Left])
			{
				ShearSpring* leftDownShear = new ShearSpring(p, p->neighbor[NodeLocation::Left], p->neighbor[NodeLocation::Down], S, R, L);
				shear_springs.push_back(leftDownShear);
			}
		}
	}

	//Add Parallel Contact Springs, avoid duplicate
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			Node* p = nodes[j*width + i];
			Node* right = p->neighbor[NodeLocation::Right];
			Node* forward = p->neighbor[NodeLocation::Down];

			if (right)
			{
				ParallelContactSpring* warpSpring = new ParallelContactSpring(p, right, Kc, R, L, Warp);
				parallel_contact_springs.push_back(warpSpring);
			}
			if (forward)
			{
				ParallelContactSpring* weftSpring = new ParallelContactSpring(p, forward, Kc, R, L, Weft);
				parallel_contact_springs.push_back(weftSpring);
			}
		}
	}

	cout << "Stretch: " << stretch_springs.size() << endl;
	cout << "Bend: " << bend_springs.size() << endl;
	cout << "Shear: " << shear_springs.size() << endl;
	cout << "Parallel: " << parallel_contact_springs.size() << endl;

	v.resize(nodes.size()*5);
	f.resize(nodes.size() * 5);

	M.resize(nodes.size() * 5, nodes.size() * 5);
	K.resize(nodes.size() * 5, nodes.size() * 5);

	v.setZero();
	f.setZero();
	M.setZero();
	K.setZero();
}

void Cloth::computeInertia()
{



}


void Cloth::computeForce()
{
	/* Initialize */
	f.setZero();
	K.setZero();
	vector<T> _K;

	for (int i = 0; i < nodes.size(); i++)
	{
		nodes[i]->compressForce = 0;
	}

	/* Compute stretching and bending, while storing compression force */
	for (int i = 0; i < stretch_springs.size(); i++)
	{
		stretch_springs[i]->solve(_K, f);
	}

	for (int i = 0; i < bend_springs.size(); i++)
	{
		bend_springs[i]->solve(_K, f);
	}

	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i]->compressForce < 0) nodes[i]->compressForce = 0;
	}

	/*
		Compute friction force here
	*/



	for (int i = 0; i < shear_springs.size(); i++)
	{
		shear_springs[i]->solve(_K, f);
	}

	for (int i = 0; i < parallel_contact_springs.size(); i++)
	{
		parallel_contact_springs[i]->solve(_K, f);
	}
}


void Cloth::draw()
{
	for (int i = 0; i < stretch_springs.size(); i++)
	{
		stretch_springs[i]->draw();
	}
}




