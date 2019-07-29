#include "cloth.h"
#include "node.h"
#include "stretch.h"
#include "bend.h"
#include "parallel_contact.h"
#include "shear.h"

#include <iostream>
#include <freeglut.h>

using namespace std;
using namespace Eigen;

Cloth::Cloth(int w, int h, double _R, double _L, double _mu, double _rho, double _Y, double _B, double _S, double _Kc, double _Kf)
{
	width = w;
	height = h;
	R = _R;
	L = _L;
	mu = _mu;
	rho = _rho;

	Y = _Y;
	B = _B;
	S = _S;
	Kc = _Kc;
	Kf = _Kf;

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

/* Fill the generalized mass matrix */
void Cloth::computeInertia()
{
	M.setZero();
	vector<T> _M;
	for (int i = 0; i < stretch_springs.size(); i++)
	{
		Node* n0 = stretch_springs[i]->node0;
		Node* n1 = stretch_springs[i]->node1;

		int index0 = n0->index * 5;
		int index1 = n1->index * 5;

		Matrix3d I = Matrix3d::Identity();

		if (stretch_springs[i]->springType == Warp)
		{
			double delta_u = n1->u - n0->u;
			Vector3d w = (n1->position - n0->position) / delta_u;

			MatrixXd m00;
			m00.resize(5, 5);
			m00.setZero();
			m00.block<3, 3>(0, 0) = 2 * I;
			m00.block<3, 1>(0, 3) = -2 * w;
			m00.block<1, 3>(3, 0) = -2 * w.transpose();
			m00(3, 3) = 2 * w.dot(w);
			m00 = rho * delta_u / 6.0*m00;
			fillBlock(_M, m00, index0, index0);

			MatrixXd m01;
			m01.resize(5, 5);
			m01.setZero();
			m01.block<3, 3>(0, 0) = I;
			m01.block<3, 1>(0, 3) = -w;
			m01.block<1, 3>(3, 0) = -w.transpose();
			m01(3, 3) = w.dot(w);
			m01 = rho * delta_u / 6.0*m01;
			fillBlock(_M, m01, index0, index1);

			MatrixXd m10;
			m10.resize(5, 5);
			m10.setZero();
			m10.block<3, 3>(0, 0) = I;
			m10.block<3, 1>(0, 3) = -w;
			m10.block<1, 3>(3, 0) = -w.transpose();
			m10(3, 3) = w.dot(w);
			m10 = rho * delta_u / 6.0*m10;
			fillBlock(_M, m10, index1, index0);

			MatrixXd m11;
			m11.resize(5, 5);
			m11.setZero();
			m11.block<3, 3>(0, 0) = 2 * I;
			m11.block<3, 1>(0, 3) = -2 * w;
			m11.block<1, 3>(3, 0) = -2 * w.transpose();
			m11(3, 3) = 2 * w.dot(w);
			m11 = rho * delta_u / 6.0*m11;
			fillBlock(_M, m11, index1, index1);
		}
		else if (stretch_springs[i]->springType == Weft)
		{
			double delta_v = n1->v - n0->v;
			Vector3d w = (n1->position - n0->position) / delta_v;

			MatrixXd m00;
			m00.resize(5, 5);
			m00.setZero();
			m00.block<3, 3>(0, 0) = 2 * I;
			m00.block<3, 1>(0, 4) = -2 * w;
			m00.block<1, 3>(4, 0) = -2 * w.transpose();
			m00(4, 4) = 2 * w.dot(w);
			m00 = rho * delta_v / 6.0*m00;
			fillBlock(_M, m00, index0, index0);

			MatrixXd m01;
			m01.resize(5, 5);
			m01.setZero();
			m01.block<3, 3>(0, 0) = I;
			m01.block<3, 1>(0, 4) = -w;
			m01.block<1, 3>(4, 0) = -w.transpose();
			m01(4, 4) = w.dot(w);
			m01 = rho * delta_v / 6.0*m01;
			fillBlock(_M, m01, index0, index1);

			MatrixXd m10;
			m10.resize(5, 5);
			m10.setZero();
			m10.block<3, 3>(0, 0) = I;
			m10.block<3, 1>(0, 4) = -w;
			m10.block<1, 3>(4, 0) = -w.transpose();
			m10(4, 4) = w.dot(w);
			m10 = rho * delta_v / 6.0*m10;
			fillBlock(_M, m10, index1, index0);

			MatrixXd m11;
			m11.resize(5, 5);
			m11.setZero();
			m11.block<3, 3>(0, 0) = 2 * I;
			m11.block<3, 1>(0, 4) = -2 * w;
			m11.block<1, 3>(4, 0) = -2 * w.transpose();
			m11(4, 4) = 2 * w.dot(w);
			m11 = rho * delta_v / 6.0*m11;
			fillBlock(_M, m11, index1, index1);
		}
	}

	M.setFromTriplets(_M.begin(), _M.end());
}

/* Fill the stiffness matrix and force vector */
void Cloth::computeForce(Vector3d gravity)
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

	for (int i = 0; i < nodes.size(); i++)
	{
		nodes[i]->getFriction(mu, Kf, _K, f);
	}

	for (int i = 0; i < stretch_springs.size(); i++)
	{
		Node* n0 = stretch_springs[i]->node0;
		Node* n1 = stretch_springs[i]->node1;

		int index0 = n0->index * 5;
		int index1 = n1->index * 5;

		if (stretch_springs[i]->springType == Warp)
		{
			double delta_u = n1->u - n0->u;
			f.segment<3>(index0) += 0.5*rho*delta_u*gravity;
			f.segment<3>(index1) += 0.5*rho*delta_u*gravity;
		}
		else if (stretch_springs[i]->springType == Weft)
		{
			double delta_v = n1->v - n0->v;
			f.segment<3>(index0) += 0.5*rho*delta_v*gravity;
			f.segment<3>(index1) += 0.5*rho*delta_v*gravity;
		}
	}

	for (int i = 0; i < shear_springs.size(); i++)
	{
		shear_springs[i]->solve(_K, f);
	}

	for (int i = 0; i < parallel_contact_springs.size(); i++)
	{
		parallel_contact_springs[i]->solve(_K, f);
	}

	K.setFromTriplets(_K.begin(), _K.end());
}


void Cloth::draw()
{
	for (int i = 0; i < stretch_springs.size(); i++)
	{
		stretch_springs[i]->draw();
	}

	for (int i = 0; i < nodes.size(); i++)
	{
		glPointSize(5.0f);
		glColor3d(0.8, 0, 0);
		glBegin(GL_POINTS);
		glVertex3d(nodes[i]->position.x(),nodes[i]->position.y(),nodes[i]->position.z());
		glEnd();
	}
}

void Cloth::step(double h)
{
	for (int i = 0; i < nodes.size(); i++)
	{
		nodes[i]->getNormal();
	}

	computeInertia();

	computeForce(Vector3d(0,0,-9.8));

	solve(h);

	/* Update */
	for (int i = 0; i < nodes.size(); i++)
	{
		if (i <width)
		{
			v.segment<3>(nodes[i]->index * 5) = Vector3d::Zero();
			v(nodes[i]->index * 5 + 3) = 0;
			v(nodes[i]->index * 5 + 4) = 0;
			continue;
		}
		
		nodes[i]->velocity = v.segment<3>(nodes[i]->index * 5);
		nodes[i]->velocityUV.x() = v(nodes[i]->index * 5 + 3);
		nodes[i]->velocityUV.y() = v(nodes[i]->index * 5 + 4);

		nodes[i]->position = nodes[i]->position + nodes[i]->velocity*h;

		nodes[i]->u = nodes[i]->u + nodes[i]->velocityUV.x()*h;
		nodes[i]->v = nodes[i]->v + nodes[i]->velocityUV.y()*h;
	}
}



void Cloth::solve(double h)
{
	SparseMatrix<double> A = M - h * h*K;
	VectorXd b = M * v + h * f;
	
	ConjugateGradient<SparseMatrix<double>, Lower | Upper> lscg;
	lscg.compute(A);
	v = lscg.solve(b);
	if (lscg.info() != Success) cout << "Failed" << endl;
}