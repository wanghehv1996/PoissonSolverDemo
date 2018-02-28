#ifndef _VISUAL_H
#define _VISUAL_H
#include <GL/glut.h>
#include "eigen3/Eigen/Dense"
#include <eigen3/Eigen/Eigen> 
#include <iostream>
#include <math.h>



class MatGraph3d{
public:
	MatGraph3d(double *mat, int dimx, int dimy, int dimz, double scale, bool renderType = true):
		_mat(mat),_dimx(dimx), _dimy(dimy),_dimz(dimz),_scale(scale),_renderType(renderType),_mask(NULL){};
	void SetMask(int *);
	void ShowResult();
private:
	int *_mask;
	double *_mat;
	int _dimx, _dimy, _dimz;
	double _scale;
	int _renderType;
};

void SetMatGraph3d(MatGraph3d *);
void InitGL(int argc, char** argv);

void display();
void keyboard(unsigned char, int, int);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int, int);
#endif
