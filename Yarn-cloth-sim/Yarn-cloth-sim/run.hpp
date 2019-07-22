#include <freeglut.h>
#include <glm/glm.hpp>
#include <Eigen/sparse>
#include <Eigen/Dense>
#include <iostream>

#include "cloth.h"
#include "camera.h"

using namespace std;
using namespace Eigen;

double mouseX = 0, mouseY = 0;
bool firstMouse = true;

Camera camera(glm::vec3(12, 11, 5), glm::vec3(0,0,1),-132,-17);
Cloth cloth(15,15,0.05,0.5);

void init()
{
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Yarn");
	glEnable(GL_DEPTH_TEST);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	
	glClearColor(0.2f, 0.2f, 0.2f, 1.0);
}

void reshape(int width, int height)
{
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0f, (GLfloat)width / (GLfloat)height, 0.05f, 100.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void display()
{
	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(camera.Position.x, camera.Position.y, camera.Position.z,
		camera.Position.x + camera.Front.x, camera.Position.y + camera.Front.y, camera.Position.z + camera.Front.z,
		camera.Up.x, camera.Up.y, camera.Up.z);

	glLineWidth(2.0);
	glColor3d(1, 0, 0);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(10, 0, 0);
	glEnd();

	glColor3d(0, 1, 0);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 10, 0);
	glEnd();

	glColor3d(0, 0, 1);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 10);
	glEnd();

	cloth.draw();

	glutPostRedisplay();
	glutSwapBuffers();
}

void mymouse(int button, int state, int x, int y)
{
	if (state = GLUT_UP)
	{
		firstMouse = true;
	}
}

void mouseMotion(int x, int y)
{
	if (firstMouse)
	{
		mouseX = x;
		mouseY = y;
		firstMouse = false;
	}

	double xoffset = x - mouseX;
	double yoffset = y - mouseY;

	mouseX = x;
	mouseY = y;

	camera.rotate(xoffset, yoffset);
	glutPostRedisplay();
}

void mykeyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'w': camera.translate(FORWARD, 0.08); break;
		case 'a': camera.translate(LEFT, 0.08); break;
		case 's': camera.translate(BACKWARD, 0.08); break;
		case 'd': camera.translate(RIGHT, 0.08); break;
	}
	glutPostRedisplay();
}

void run(int argc, char* argv[])
{
	glutInit(&argc, (char**)argv);
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(mykeyboard);
	glutMouseFunc(mymouse);
	glutMotionFunc(mouseMotion);
	glutMainLoop();
}