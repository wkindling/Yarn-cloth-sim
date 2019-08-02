#include "cloth.h"
#include "node.h"
#include "stretch.h"
#include "bend.h"
#include "parallel_contact.h"
#include "shear.h"

#include "Mosek/QuadProgMosek.h"

#include <iostream>
#include <freeglut.h>

#define MOSEK_SOLVE

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

	DoF = 0;

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

			node->node_index = nodes.size(); //Record index for later global matrix
			nodes.push_back(node);
		}
	}

	int cross_count = 0;
	for (int j = 0; j < height; j++) //Record the neighbors of each node
	{
		for (int i = 0; i < width; i++)
		{
			Node* node = nodes[j*width + i];
			node->neighbor[NodeLocation::Left] = i >= 1 ? nodes[j*width + i - 1] : NULL;
			node->neighbor[NodeLocation::Right] = i < width - 1 ? nodes[j*width + i + 1] : NULL;
			node->neighbor[NodeLocation::Up] = j >= 1 ? nodes[(j - 1)*width + i] : NULL;
			node->neighbor[NodeLocation::Down] = j < height - 1 ? nodes[(j + 1)*width + i] : NULL;

			for (int i = 0; i < 4; i++)
			{
				if (node->neighbor[i] == NULL) node->onBorder = true;	
			}

			if (!node->onBorder)
			{
				node->cross_index = cross_count;
				cross_count++;
			}
		}
	}

	/* Apply fixed constraint */
	for (int i = 0; i < nodes.size(); i++)
	{
		if (i == width - 1 || i == 0)
		{
			nodes[i]->fixed = true;
		}
	}

	DoF = 3 * nodes.size() + 2 * cross_count;

	cout << "DOF : " << DoF << endl;

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

	v.resize(DoF);
	f.resize(DoF);

	M.resize(DoF, DoF);
	K.resize(DoF, DoF);

	v.setZero();
	f.setZero();
	M.setZero();
	K.setZero();
}

/* Fill the generalized mass matrix */
void Cloth::computeInertia(vector<T>& _K)
{
	M.setZero();
	vector<T> _M;

	for (int i = 0; i < stretch_springs.size(); i++)
	{
		Node* n0 = stretch_springs[i]->node0;
		Node* n1 = stretch_springs[i]->node1;

		int node_index0 = n0->node_index * 3;
		int node_index1 = n1->node_index * 3;

		Matrix3d I = Matrix3d::Identity();

		if (stretch_springs[i]->springType == Warp)
		{
			double delta_u = abs(n1->u - n0->u);
			double coeff = 1.0 / 6.0*rho*delta_u;

			Vector3d w = (n1->position - n0->position) / delta_u;

			/* Lagrange - Lagrange Part */
			MatrixXd L0L0 = coeff * 2 * I;
			MatrixXd L0L1 = coeff * I;
			MatrixXd L1L0 = coeff * I;
			MatrixXd L1L1 = coeff * 2 * I;

			fillGlobalInertia(_M, L0L0, node_index0, node_index0);
			fillGlobalInertia(_M, L0L1, node_index0, node_index1);
			fillGlobalInertia(_M, L1L0, node_index1, node_index0);
			fillGlobalInertia(_M, L1L1, node_index1, node_index1);

			fillGlobalInertia(_K, L0L0, node_index0, node_index0);
			fillGlobalInertia(_K, L0L1, node_index0, node_index1);
			fillGlobalInertia(_K, L1L0, node_index1, node_index0);
			fillGlobalInertia(_K, L1L1, node_index1, node_index1);


			/* Euler Part */
			if (!n0->onBorder)
			{
				int cross_index0 = nodes.size() * 3 + n0->cross_index * 2;

				MatrixXd L0E0 = -2 * coeff* w;
				MatrixXd E0L0 = -2 * coeff* w.transpose();
				fillGlobalInertia(_M, L0E0, node_index0, cross_index0);
				fillGlobalInertia(_M, E0L0, cross_index0, node_index0);

				fillGlobalInertia(_K, L0E0, node_index0, cross_index0);
				fillGlobalInertia(_K, E0L0, cross_index0, node_index0);

				_M.push_back(T(cross_index0, cross_index0, 2 * coeff*w.dot(w)));
				_K.push_back(T(cross_index0, cross_index0, 2 * coeff*w.dot(w)));
			}

			if (!n0->onBorder && !n1->onBorder)
			{
				int cross_index0 = nodes.size() * 3 + n0->cross_index * 2;
				int cross_index1 = nodes.size() * 3 + n1->cross_index * 2;

				MatrixXd L0E1 = -coeff * w;
				MatrixXd E0L1 = -coeff * w.transpose();
				fillGlobalInertia(_M, L0E1, node_index0, cross_index1);
				fillGlobalInertia(_M, E0L1, cross_index0, node_index1);

				fillGlobalInertia(_K, L0E1, node_index0, cross_index1);
				fillGlobalInertia(_K, E0L1, cross_index0, node_index1);

				_M.push_back(T(cross_index0, cross_index1, coeff*w.dot(w)));
				_K.push_back(T(cross_index0, cross_index1, coeff*w.dot(w)));

				MatrixXd L1E0 = -coeff * w;
				MatrixXd E1L0 = -coeff * w.transpose();
				fillGlobalInertia(_M, L1E0, node_index1, cross_index0);
				fillGlobalInertia(_M, E1L0, cross_index1, node_index0);

				fillGlobalInertia(_K, L1E0, node_index1, cross_index0);
				fillGlobalInertia(_K, E1L0, cross_index1, node_index0);

				_M.push_back(T(cross_index1, cross_index0, coeff*w.dot(w)));
				_K.push_back(T(cross_index1, cross_index0, coeff*w.dot(w)));
			}

			if (!n1->onBorder)
			{
				int cross_index1 = nodes.size() * 3 + n1->cross_index * 2;

				MatrixXd L1E1 = -2 * coeff*w;
				MatrixXd E1L1 = -2 * coeff*w.transpose();
				fillGlobalInertia(_M, L1E1, node_index1, cross_index1);
				fillGlobalInertia(_M, E1L1, cross_index1, node_index1);

				fillGlobalInertia(_K, L1E1, node_index1, cross_index1);
				fillGlobalInertia(_K, E1L1, cross_index1, node_index1);

				_M.push_back(T(cross_index1, cross_index1, 2 * coeff*w.dot(w)));
				_K.push_back(T(cross_index1, cross_index1, 2 * coeff*w.dot(w)));
			}
		}
		else if (stretch_springs[i]->springType == Weft)
		{
			double delta_v = abs(n1->v - n0->v);
			double coeff = 1.0 / 6.0*rho*delta_v;

			Vector3d w = (n1->position - n0->position) / delta_v;

			/* Lagrange-Lagrange Part */
			MatrixXd L0L0 = coeff * 2 * I;
			MatrixXd L0L1 = coeff * I;
			MatrixXd L1L0 = coeff * I;
			MatrixXd L1L1 = coeff * 2 * I;

			fillGlobalInertia(_M, L0L0, node_index0, node_index0);
			fillGlobalInertia(_M, L0L1, node_index0, node_index1);
			fillGlobalInertia(_M, L1L0, node_index1, node_index0);
			fillGlobalInertia(_M, L1L1, node_index1, node_index1);

			fillGlobalInertia(_K, L0L0, node_index0, node_index0);
			fillGlobalInertia(_K, L0L1, node_index0, node_index1);
			fillGlobalInertia(_K, L1L0, node_index1, node_index0);
			fillGlobalInertia(_K, L1L1, node_index1, node_index1);

			/* Euler Part */
			if (!n0->onBorder)
			{
				int cross_index0 = nodes.size() * 3 + n0->cross_index * 2;

				MatrixXd L0E0 = -2 * coeff*w;
				MatrixXd E0L0 = -2 * coeff*w.transpose();
				fillGlobalInertia(_M, L0E0, node_index0, cross_index0 + 1);
				fillGlobalInertia(_M, E0L0, cross_index0 + 1, node_index1);

				fillGlobalInertia(_K, L0E0, node_index0, cross_index0 + 1);
				fillGlobalInertia(_K, E0L0, cross_index0 + 1, node_index1);

				_M.push_back(T(cross_index0 + 1, cross_index0 + 1, 2 * coeff*w.dot(w)));
				_K.push_back(T(cross_index0 + 1, cross_index0 + 1, 2 * coeff*w.dot(w)));
			}

			if (!n0->onBorder && !n1->onBorder)
			{
				int cross_index0 = nodes.size() * 3 + n0->cross_index * 2;
				int cross_index1 = nodes.size() * 3 + n1->cross_index * 2;

				MatrixXd L0E1 = -coeff * w;
				MatrixXd E0L1 = -coeff * w.transpose();
				fillGlobalInertia(_M, L0E1, node_index0, cross_index1 + 1);
				fillGlobalInertia(_M, E0L1, cross_index0 + 1, node_index1);

				fillGlobalInertia(_K, L0E1, node_index0, cross_index1 + 1);
				fillGlobalInertia(_K, E0L1, cross_index0 + 1, node_index1);

				_M.push_back(T(cross_index0 + 1, cross_index1 + 1, coeff*w.dot(w)));
				_K.push_back(T(cross_index0 + 1, cross_index1 + 1, coeff*w.dot(w)));

				MatrixXd L1E0 = -coeff * w;
				MatrixXd E1L0 = -coeff * w.transpose();
				fillGlobalInertia(_M, L1E0, node_index1, cross_index0 + 1);
				fillGlobalInertia(_M, E1L0, cross_index1 + 1, node_index0);

				fillGlobalInertia(_K, L1E0, node_index1, cross_index0 + 1);
				fillGlobalInertia(_K, E1L0, cross_index1 + 1, node_index0);

				_M.push_back(T(cross_index1 + 1, cross_index0 + 1, coeff*w.dot(w)));
				_K.push_back(T(cross_index1 + 1, cross_index0 + 1, coeff*w.dot(w)));
			}

			if (!n1->onBorder)
			{
				int cross_index1 = nodes.size() * 3 + n1->cross_index * 2;

				MatrixXd L1E1 = -2 * coeff*w;
				MatrixXd E1L1 = -2 * coeff*w.transpose();
				fillGlobalInertia(_M, L1E1, node_index1, cross_index1 + 1);
				fillGlobalInertia(_M, E1L1, cross_index1 + 1, node_index1);

				fillGlobalInertia(_K, L1E1, node_index1, cross_index1 + 1);
				fillGlobalInertia(_K, E1L1, cross_index1 + 1, node_index1);

				_M.push_back(T(cross_index1 + 1, cross_index1 + 1, 2 * coeff*w.dot(w)));
				_K.push_back(T(cross_index1 + 1, cross_index1 + 1, 2 * coeff*w.dot(w)));
			}
		}
	}
	M.setFromTriplets(_M.begin(), _M.end());
}

/* Fill the stiffness matrix and force vector */
void Cloth::computeForce(Vector3d gravity, double h)
{
	/* Initialize */
	f.setZero();
	K.setZero();
	vector<T> _K;

	computeInertia(_K);

	for (int i = 0; i < nodes.size(); i++)
	{
		nodes[i]->compressForce = 0;
	}

	/* Compute stretching and bending, while storing compression force */
	for (int i = 0; i < stretch_springs.size(); i++)
	{
		stretch_springs[i]->solve(_K, f, nodes.size(), h);
	}

	for (int i = 0; i < bend_springs.size(); i++)
	{
		bend_springs[i]->solve(_K, f, nodes.size(), h);
	}

	/* Get Friction force */
	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i]->compressForce < 0) nodes[i]->compressForce = 0;
	}

	for (int i = 0; i < nodes.size(); i++)
	{
		nodes[i]->getFriction(mu, Kf, _K, f, nodes.size(), h);
	}

	/* Gravity Force */
	for (int i = 0; i < stretch_springs.size(); i++)
	{
		Node* n0 = stretch_springs[i]->node0;
		Node* n1 = stretch_springs[i]->node1;

		int node_index0 = n0->node_index * 3;
		int node_index1 = n1->node_index * 3;

		if (stretch_springs[i]->springType == Warp)
		{
			double delta_u = n1->u - n0->u;
			f.segment<3>(node_index0) += 0.5*rho*delta_u*gravity;
			f.segment<3>(node_index1) += 0.5*rho*delta_u*gravity;

			if (!n0->onBorder)
			{
				int cross_index0 = nodes.size() * 3 + n0->cross_index * 2;
				double Fu0 = -rho * gravity.dot((n0->position + n1->position) / 2.0);
				f(cross_index0) += Fu0;

				MatrixXd Fu0dx0 = -0.5*rho*gravity.transpose();
				MatrixXd Fu0dx1 = -0.5*rho*gravity.transpose();
				fillGlobalStiffness(_K, Fu0dx0, cross_index0, node_index0, h);
				fillGlobalStiffness(_K, Fu0dx1, cross_index0, node_index1, h);

				MatrixXd Fx0du0 = -0.5*rho*gravity;
				MatrixXd Fx1du0 = -0.5*rho*gravity;
				fillGlobalStiffness(_K, Fx0du0, node_index0, cross_index0, h);
				fillGlobalStiffness(_K, Fx1du0, node_index1, cross_index0, h);
			}

			if (!n1->onBorder)
			{
				int cross_index1 = nodes.size() * 3 + n1->cross_index * 2;
				double Fu1 = rho * gravity.dot((n0->position + n1->position) / 2.0);
				f(cross_index1) += Fu1;

				MatrixXd Fu1dx0 = 0.5*rho*gravity.transpose();
				MatrixXd Fu1dx1 = 0.5*rho*gravity.transpose();
				fillGlobalStiffness(_K, Fu1dx0, cross_index1, node_index0, h);
				fillGlobalStiffness(_K, Fu1dx1, cross_index1, node_index1, h);

				MatrixXd Fx0du1 = 0.5*rho*gravity;
				MatrixXd Fx1du1 = 0.5*rho*gravity;
				fillGlobalStiffness(_K, Fx0du1, node_index0, cross_index1, h);
				fillGlobalStiffness(_K, Fx1du1, node_index1, cross_index1, h);
			}
		}
		else if (stretch_springs[i]->springType == Weft)
		{
			double delta_v = n1->v - n0->v;
			f.segment<3>(node_index0) += 0.5*rho*delta_v*gravity;
			f.segment<3>(node_index1) += 0.5*rho*delta_v*gravity;

			if (!n0->onBorder)
			{
				int cross_index0 = nodes.size() * 3 + n0->cross_index * 2;
				double Fv0 = -rho * gravity.dot((n0->position + n1->position) / 2.0);
				f(cross_index0 + 1) += Fv0;

				MatrixXd Fv0dx0 = -0.5*rho*gravity.transpose();
				MatrixXd Fv0dx1 = -0.5*rho*gravity.transpose();
				fillGlobalStiffness(_K, Fv0dx0, cross_index0 + 1, node_index0, h);
				fillGlobalStiffness(_K, Fv0dx1, cross_index0 + 1, node_index1, h);

				MatrixXd Fx0dv0 = -0.5*rho*gravity;
				MatrixXd Fx1dv0 = -0.5*rho*gravity;
				fillGlobalStiffness(_K, Fx0dv0, node_index0, cross_index0 + 1, h);
				fillGlobalStiffness(_K, Fx1dv0, node_index1, cross_index0 + 1, h);
			}

			if (!n1->onBorder)
			{
				int cross_index1 = nodes.size() * 3 + n1->cross_index * 2;
				double Fv1 = rho * gravity.dot((n0->position + n1->position) / 2.0);
				f(cross_index1 + 1) += Fv1;

				MatrixXd Fv1dx0 = 0.5*rho*gravity.transpose();
				MatrixXd Fv1dx1 = 0.5*rho*gravity.transpose();
				fillGlobalStiffness(_K, Fv1dx0, cross_index1 + 1, node_index0, h);
				fillGlobalStiffness(_K, Fv1dx1, cross_index1 + 1, node_index1, h);

				MatrixXd Fx0dv1 = 0.5*rho*gravity;
				MatrixXd Fx1dv1 = 0.5*rho*gravity;
				fillGlobalStiffness(_K, Fx0dv1, node_index0, cross_index1 + 1, h);
				fillGlobalStiffness(_K, Fx1dv1, node_index1, cross_index1 + 1, h);
			}
		}
	}

	for (int i = 0; i < shear_springs.size(); i++)
	{
		shear_springs[i]->solve(_K, f, h);
	}

	for (int i = 0; i < parallel_contact_springs.size(); i++)
	{
		parallel_contact_springs[i]->solve(_K, f, nodes.size(), h);
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
		glVertex3d(nodes[i]->position.x()*1e3,nodes[i]->position.y()*1e3,nodes[i]->position.z()*1e3);
		glEnd();
	}
}

void Cloth::step(double h)
{
	for (int i = 0; i < nodes.size(); i++)
	{
		nodes[i]->getNormal();
	}

	computeForce(Vector3d(0, 0, -9.8), h);
	
	applyConstraint();

	solve(h);

	/* Update */
	for (int i = 0; i < nodes.size(); i++)
	{
		nodes[i]->velocity = v.segment<3>(nodes[i]->node_index * 3);
		nodes[i]->position = nodes[i]->position + nodes[i]->velocity*h;

		if (!nodes[i]->onBorder)
		{
			nodes[i]->velocityUV.x() = v(nodes.size() * 3 + nodes[i]->cross_index * 2);
			nodes[i]->velocityUV.y() = v(nodes.size() * 3 + nodes[i]->cross_index * 2 + 1);

			nodes[i]->u = nodes[i]->u + nodes[i]->velocityUV.x()*h;
			nodes[i]->v = nodes[i]->v + nodes[i]->velocityUV.y()*h;
		}
	}
}

void Cloth::solve(double h)
{
	SparseMatrix<double> A = K;
	VectorXd b = M * v + h * f;
	
#ifdef EIGEN_SOLVE

	LeastSquaresConjugateGradient<SparseMatrix<double>> lscg;
	lscg.compute(A);
	v = lscg.solve(b);

	if (lscg.info() != Success)
	{;
		cout << "Failed " << lscg.info() << endl;
	}

#endif

#ifdef MOSEK_SOLVE
	
	bool success = mosekSolve(A, -b, Aeq, beq, Aineq, bineq, v);

#endif
}

bool Cloth::mosekSolve(const SparseMatrix<double>& MDK, const VectorXd& b,
	const SparseMatrix<double>& Aeq, const VectorXd& beq,
	const SparseMatrix<double>& Aineq, const VectorXd& bineq,
	VectorXd& v)
{
	QuadProgMosek *program = new QuadProgMosek();
	double inf = numeric_limits<double>::infinity();

	VectorXd xl;
	VectorXd xu;
	xl.setConstant(b.size(), -inf);
	xu.setConstant(b.size(), inf);

	program->setNumberOfVariables(b.size());
	program->setNumberOfEqualities(beq.size());
	program->setNumberOfInequalities(bineq.size());

	program->setObjectiveMatrix(MDK);
	program->setObjectiveVector(b);

	program->setInequalityMatrix(Aineq);
	program->setInequalityVector(bineq);

	program->setEqualityMatrix(Aeq);
	program->setEqualityVector(beq);

	bool success = program->solve();

	v = program->getPrimalSolution();

	return success;
}

void Cloth::applyConstraint()
{
	vector<T> _Aeq;
	vector<T> _Aineq;
	vector<pair<int, double>> _beq;
	vector<pair<int, double>> _bineq;
	
	int eqsize = 0;
	int ineqsize = 0;

	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i]->fixed)
		{
			_Aeq.push_back(T(eqsize, nodes[i]->node_index * 3, 1));
			_beq.push_back(make_pair(eqsize, nodes[i]->velocity.x()));
			eqsize++;

			_Aeq.push_back(T(eqsize, nodes[i]->node_index * 3 + 1, 1));
			_beq.push_back(make_pair(eqsize, nodes[i]->velocity.y()));
			eqsize++;
			
			_Aeq.push_back(T(eqsize, nodes[i]->node_index * 3 + 2, 1));
			_beq.push_back(make_pair(eqsize, nodes[i]->velocity.z()));
			eqsize++;

			if (!nodes[i]->onBorder)
			{
				_Aeq.push_back(T(eqsize, nodes[i]->cross_index * 2 + nodes.size() * 3, 1));
				_beq.push_back(make_pair(eqsize, nodes[i]->velocityUV.x()));
				eqsize++;

				_Aeq.push_back(T(eqsize, nodes[i]->cross_index * 2 + nodes.size() * 3 + 1, 1));
				_beq.push_back(make_pair(eqsize, nodes[i]->velocityUV.y()));
				eqsize++;
			}
		}
	}
	
	Aeq.resize(eqsize, DoF);
	Aineq.resize(ineqsize, DoF);
	beq.resize(eqsize);
	bineq.resize(ineqsize);

	Aeq.setFromTriplets(_Aeq.begin(), _Aeq.end());
	Aineq.setFromTriplets(_Aineq.begin(), _Aineq.end());

	beq.setZero();
	bineq.setZero();

	for (int i = 0; i < _beq.size(); i++)
	{
		beq(_beq[i].first) = _beq[i].second;
	}

	for (int i = 0; i < _bineq.size(); i++)
	{
		bineq(_bineq[i].first) = _bineq[i].second;
	}
}