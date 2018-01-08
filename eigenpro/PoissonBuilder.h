#ifndef _POISSON_BUILDER_H
#define _POISSON_BUILDER_H
#include "eigen3/Eigen/Dense"
#include <eigen3/Eigen/Eigen> 
#include <eigen3/Eigen/SparseCore>
#include <vector>
#include <stdio.h>

#define PB_BOUND_Dirichlet	1
#define PB_BOUND_Neumann	2

#define PB_BOUND_LEFT	4
#define PB_BOUND_RIGHT	8
#define PB_BOUND_UP		16
#define PB_BOUND_DOWN	32



class PoissonBuilder{
public:
	//Initialize
	PoissonBuilder(int dimx, int dimy, double dx):
		_dimx(dimx), _dimy(dimy), _dx(dx){};
	
	//Set ForceFunction
	int SetForce(double (*forceFunc)(int, int, int, int));

	//Set Boundary Condition & Boundary Function
	int SetBoundary(int boundarytype, double (*boundaryFunc)(int,int));

	//Get the Laplacian Matrix(A)
	Eigen::SparseMatrix<double> GetLaplacianMatrix();

	//Get the Force Vector(b)
	Eigen::VectorXd GetForceVector();

	void SetDebugInfo(bool debugInfo){_debugInfo = debugInfo;};

private:
	//dimension of x, y
	int _dimx, _dimy;

	//delta x
	double _dx;

	//force function
	double (*_forceFunc)(int, int, int, int);//the forcefunction

	//type of each boundary(PB_BOUND_Dirichlet / PB_BOUND_Neumann)
	int _boundaryType[4] = {0,0,0,0};

	//function of each boundary
	double (*_boundaryFunc[4])(int, int) = {NULL, NULL, NULL, NULL};

	//open debug infomation or not
	bool _debugInfo = true;

};


double zeroff(int, int, int, int);
double sinpowff(int, int, int, int);
double sinff(int, int, int, int);

double sinpowbc(int i, int total);
double sinbc(int i, int total);

#endif