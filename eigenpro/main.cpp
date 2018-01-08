#include "eigen3/Eigen/Dense"
#include <eigen3/Eigen/Eigen> 
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <iostream>
#include <stdlib.h>
using namespace std;
using namespace Eigen;

#include "visual.h"
#include "PoissonBuilder.h"

/*
------> x
|
|
v y

Matrix(y,x)
y = row num
x = col num
*/

void GetDivergence(MatrixXd &result, MatrixXd originMapX, MatrixXd originMapY, double invdx, int boundStrategy){
	result.resize(originMapX.rows(), originMapX.cols());
	for(int i=0;i<originMapX.rows();i++)
		for(int j=0;j<originMapX.cols();j++)
		{
			double yplus = (i<originMapY.rows()-1)?
				originMapY(i + 1, j):(boundStrategy?originMapY(i,j):0);
			double yminus = (i>0)?
				originMapY(i - 1, j):(boundStrategy?originMapY(i,j):0);

			double xplus = (j<originMapX.cols()-1)?
				originMapX(i, j+1):(boundStrategy?originMapX(i,j):0);
			double xminus = (j>0)?
				originMapX(i,j - 1):(boundStrategy?originMapX(i,j):0);

			
			double dY = yplus-yminus;
			double dX = xplus-xminus;
			result(i,j) = 0.5*(yplus-yminus)*invdx + 0.5*(xplus-xminus)*invdx;
		}
	//return result;
}

//
void GetGradient(MatrixXd &resultX, MatrixXd &resultY, MatrixXd originMap, double invdx, int boundStrategy){
	resultX.resize(originMap.rows(), originMap.cols());
	resultY.resize(originMap.rows(), originMap.cols());
	for(int i=0;i<originMap.rows();i++)
		for(int j=0;j<originMap.cols();j++)
		{
			double yplus = (i<originMap.rows()-1)?
				originMap(i + 1, j):(boundStrategy?originMap(i,j):0);
			double yminus = (i>0)?
				originMap(i - 1, j):(boundStrategy?originMap(i,j):0);

			double xplus = (j<originMap.cols()-1)?
				originMap(i, j+1):(boundStrategy?originMap(i,j):0);
			double xminus = (j>0)?
				originMap(i,j - 1):(boundStrategy?originMap(i,j):0);

			
			double dY = yplus-yminus;
			double dX = xplus-xminus;
			resultY(i, j) = 0.5*(yplus-yminus)*invdx;
			resultX(i, j) = 0.5*(xplus-xminus)*invdx;
		}
	//return result;
}

//GetLaplacian
//
void GetLaplacian(MatrixXd &result, MatrixXd originMap, double invdx, int boundStrategy)
{	
	result.resize(originMap.rows(), originMap.cols());
	for(int i=0;i<originMap.rows();i++)
		for(int j=0;j<originMap.cols();j++)
		{
			double yplus = (i<originMap.rows()-1)?
				originMap(i + 1, j):(boundStrategy?originMap(i - 1,j):0);
			double yminus = (i>0)?
				originMap(i - 1, j):(boundStrategy?originMap(i + 1,j):0);

			double xplus = (j<originMap.cols()-1)?
				originMap(i, j + 1):(boundStrategy?originMap(i,j - 1):0);
			double xminus = (j>0)?
				originMap(i, j - 1):(boundStrategy?originMap(i,j + 1):0);
			
			//with half d
			double ddY = yplus+yminus-originMap(i,j)*2;
			double ddX = xplus+xminus-originMap(i,j)*2;
			result(i,j) = (ddY+ddX)*invdx*invdx;
		}
}


MatrixXd solveLalacianProblem(){
	int xs = 50;
	int ys = 50;
	PoissonBuilder poissonBuilder(xs, ys, 1);
	
	poissonBuilder.SetForce(sinpowff);
	poissonBuilder.SetBoundary(PB_BOUND_Dirichlet | 
		PB_BOUND_UP | PB_BOUND_LEFT |
		PB_BOUND_DOWN | PB_BOUND_RIGHT,
		sinpowbc);
	
	//solve Ax=b
	//you can use different solver here
	SparseMatrix<double> A = poissonBuilder.GetLaplacianMatrix();
	//cout<<endl<<"A1 = \n"<<A<<endl;
	Eigen::SimplicialCholesky<SparseMatrix<double> > chol(A);

	VectorXd bv = poissonBuilder.GetForceVector();

	Eigen::MatrixXd bm(ys,xs);
	for(int i=0;i<ys;i++)
		for(int j=0;j<xs;j++){
			bm(i,j) = bv(i*xs+j);
		}
	//cout<<endl<< "b1 = \n"<<bm<<endl;

	Eigen::VectorXd x = chol.solve(bv);
	MatrixXd xm(ys,xs);
	for(int i=0;i<ys;i++)
		for(int j=0;j<xs;j++){
			xm(i,j) = x(i*xs+j);
		}
	return xm;
}

int main(int argc, char** argv){

	//solve Ax = b; 
	//arg:row num, col num, forcefunction(b)
	MatrixXd result = solveLalacianProblem();
	//arg: matrix, horizontal scale, vertical scale, mesh or line?
	MatGraph * mat = new MatGraph(result, 0.2,0.01,true);
	SetMatGraph(mat);
	InitGL(argc,argv);
	return 0;
}
