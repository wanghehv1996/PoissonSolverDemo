#ifndef _POISSON_BUILDER_H
#define _POISSON_BUILDER_H
#include "eigen3/Eigen/Dense"
#include <eigen3/Eigen/Eigen> 
#include <eigen3/Eigen/SparseCore>
#include <vector>
#include <stdio.h>

//sqrt(Pi)
#define SQRTPI 1.77245385091

#define PB_BOUND_Dirichlet	1
#define PB_BOUND_Neumann	2

#define PB_BOUND_X0	4
#define PB_BOUND_X1	8
#define PB_BOUND_Y0	16
#define PB_BOUND_Y1	32



class PoissonBuilder{
public:
	//Initialize
	PoissonBuilder(int dimy, int dimx, double dx):
		_dimy(dimy), _dimx(dimx), _dx(dx)
		{ 
			int totalDim = dimx*dimy;
			_holeMaskMap = new int[totalDim]; //all zero
			_holeBoundaryMap = new double[totalDim]; //all zero
		};

	~PoissonBuilder(){
		delete _holeMaskMap;
		delete _holeBoundaryMap;
	}


	
	//Set ForceFunction
	int SetForce(double (*forceFunc)(int, int, int, int));

	//Set Boundary Condition & Boundary Function
	int SetBoundary(int boundarytype, double (*boundaryFunc)(int,int,int,int));

	int SetHole(int y, int x, int type, double v);

	int SetGradientX(double (*forceFunc)(int, int, int, int));

	int SetGradientY(double (*forceFunc)(int, int, int, int));

	//Get the Laplacian Matrix(A)
	Eigen::SparseMatrix<double> GetLaplacianMatrix();

	//Get the Force Vector(b)
	Eigen::VectorXd GetForceVector();

	int* GetHoleMask();

	void SetDebugInfo(bool debugInfo){_debugInfo = debugInfo;};

private:
	//dimension of x, y
	int _dimx, _dimy;

	//delta x
	double _dx;

	//force function
	double (*_forceFunc)(int, int, int, int);//the forcefunction

	//gradient function
	double (*_gradientXFunc)(int, int, int, int);//the forcefunction
	double (*_gradientYFunc)(int, int, int, int);//the forcefunction

	//type of each boundary(PB_BOUND_Dirichlet / PB_BOUND_Neumann)
	int _boundaryType[4] = {0,0,0,0};

	//function of each boundary
	double (*_boundaryFunc[4])(int, int, int, int) = {NULL, NULL, NULL, NULL};

	//open debug infomation or not
	bool _debugInfo = true;

	//0 = normal grid; 1 = Dirichlet; 2 = Neumann;
	int *_holeMaskMap = NULL;

	//if mask = 1/2, save exact value here
	//if mask = 0 && beside a Neumann boundary, save gradient value here
	//if mask = 0 && !beside a Neumann boundary, save nothing here
	double *_holeBoundaryMap = NULL;
};

double sinpowof(int yi, int xi, int dimy, int dimx);
double sinpowff(int, int, int, int);

double sinpowbc(int, int, int, int);

double sinpowgx(int yi, int xi, int dimy, int dimx);
double sinpowgy(int yi, int xi, int dimy, int dimx);


//z = x*x-y*y
double plusff(int, int, int, int);

double plusbc(int, int, int, int);

double plusgx(int yi, int xi, int dimy, int dimx);
double plusgy(int yi, int xi, int dimy, int dimx);
#endif