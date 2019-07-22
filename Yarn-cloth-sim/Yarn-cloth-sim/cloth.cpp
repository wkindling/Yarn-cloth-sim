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

	cout << width << " " << height << " " << R << " " << L << endl ;
	for (int j = 0; j < height; j++) //u for warp (x)
	{ 
		for (int i = 0; i < width; i++) // v for weft (y)
		{
			Node* node = new Node(Vector3d(i*L,j*L,0),i*L,j*L);

			cout << Vector3d(i*L, j*L, 0).transpose() << endl;

			if ((i + j) % 2 == 0)
			{
				node->whichUp = Warp;
			}
			else node->whichUp = Weft;

			/*
			if (i == 0 || i == width - 1 || j == 0 || j == height - 1)
			{
				node->whichUp = Normal;
				node->normal_index = nindex++;
				node->cross_index = -1;
			}
			else
			{
				node->normal_index = -1;
				node->cross_index = cindex++;
			}
			*/

			nodes.push_back(node);
		}
	}

	/*TODO: Add Different Springs */
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			Node* p = nodes[j*width + i];
			Node* right = (i < width - 1) ? nodes[j*width + i + 1] : NULL;
			Node* forward = (j < height - 1) ? nodes[(j + 1)*width + i] : NULL;
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





}

Cloth::~Cloth() {}



void Cloth::draw()
{
	for (int i = 0; i < stretch_springs.size(); i++)
	{
		stretch_springs[i]->draw();
	}
}

