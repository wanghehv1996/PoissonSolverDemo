#ifndef _VISUAL_H
#define _VISUAL_H
#include <GL/glut.h>
#include "eigen3/Eigen/Dense"
#include <eigen3/Eigen/Eigen> 
#include <iostream>
#include <math.h>



class MatGraph{
public:
	MatGraph(Eigen::MatrixXd mat, double horizontalScale, double verticalScale, bool renderType = true):
		_mat(mat),_horizontalScale(horizontalScale),_verticalScale(verticalScale),_renderType(renderType){};
	void ShowResult();
private:
	Eigen::MatrixXd _mat;
	double _horizontalScale;
	double _verticalScale;
	int _renderType;
};

void SetMatGraph(MatGraph *);
void InitGL(int argc, char** argv);

void display();
void keyboard(unsigned char, int, int);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int, int);
#endif
