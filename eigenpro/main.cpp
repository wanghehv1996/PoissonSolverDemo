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

#include "visual3d.h"
#include "PoissonBuilder.h"

double absoluteError(double *mat1, double *mat2, int size){
	double totalerror = 0;
	for(int i=0;i<size;i++){
		totalerror += (mat1[i]-mat2[i]) * (mat1[i]-mat2[i]);
	}
	return totalerror;
}

int main(int argc, char** argv){

	//1. Build the problem
	int dimx = 10;
	int dimy = 10;
	int dimz = 10;
	int gridNum = dimx*dimy*dimz;

	PoissonBuilder poissonBuilder(dimz, dimy, dimx, 1);
	poissonBuilder.SetForce(plusff);
	poissonBuilder.SetGradientX(plusgx);
	poissonBuilder.SetGradientY(plusgy);
	poissonBuilder.SetGradientZ(plusgz);
	poissonBuilder.SetBoundary(PB_BOUND_Dirichlet | 
		PB_BOUND_X0 | PB_BOUND_Y0 | PB_BOUND_Z0 |
		PB_BOUND_X1 | PB_BOUND_Y1 | PB_BOUND_Z1,
		plusbc);
	// poissonBuilder.SetBoundary(PB_BOUND_Neumann | 
	// 	PB_BOUND_X0 | PB_BOUND_X1,
	// 	plusgx);	
	// poissonBuilder.SetBoundary(PB_BOUND_Neumann | 
	// 	PB_BOUND_Y0 | PB_BOUND_Y1,
	// 	plusgy);
	// poissonBuilder.SetBoundary(PB_BOUND_Neumann | 
	// 	PB_BOUND_Z0 | PB_BOUND_Z1,
	// 	plusgz);

	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -SQRTPI,SQRTPI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -SQRTPI,SQRTPI);

	//2. Set Some Holes
	for(int I=0;I<1;I++)
	for(int J=0;J<1;J++)
	for(int K=0;K<1;K++)
		for(int i=2;i<8;i++)
		for(int j=2;j<8;j++)
		for(int k=2;k<8;k++){
			poissonBuilder.SetHole(I*10+i,J*10+j,K*10+k,PB_BOUND_Neumann,i+j+k);
			gridNum--;
		}

	// for(int i = 3;i<7;i++)
	// 	for(int j = 3;j<7;j++)
	// 		for(int k = 3;k<7;k++){
	// 			poissonBuilder.SetHole(i,j,k,2,i+j+k);
	// 			gridNum--;
	// 		}

	//3. Solve the Problem Ax=b
	//3.1. Get Laplacian Matrix A
	printf("GetLaplacianMatrix\n");
	SparseMatrix<double> A = poissonBuilder.GetLaplacianMatrix();
	//cout<<endl<<"A1 = \n"<<A<<endl;

	//3.2. Try Different Solver
	Eigen::SparseLU<SparseMatrix<double> > chol(A);
	/*
	These solvers cannot solve the problem correctly:
	1. Eigen::SimplicialCholesky<SparseMatrix<double> > chol(A);
	2. Eigen::SimplicialLLT<SparseMatrix<double> > chol(A);
	3. Eigen::SimplicialLDLT<SparseMatrix<double> > chol(A);
	4. Eigen::ConjugateGradient<SparseMatrix<double> > chol(A);
	
	These solvers can solve the problem correctly
	1. SparseQR-succeed, slow
		Eigen::SparseQR<SparseMatrix<double> , AMDOrdering < int > > chol;
		chol.compute(A);
	2. SparseLU-succeed, quick
		Eigen::SparseLU<SparseMatrix<double> > chol(A);
	3. LeastSquaresConjugateGradient-succeed, quick
		Eigen::LeastSquaresConjugateGradient <SparseMatrix<double> > chol(A);
	4. BiCGSTAB-succeed, quick
		Eigen::BiCGSTAB<SparseMatrix<double> > chol(A);
	*/

	//3.3. Get Force Vector
	printf("GetForceVector\n");
	VectorXd bv = poissonBuilder.GetForceVector();
	double bm[dimx*dimy*dimz];
	for(int i=0;i<dimx*dimy*dimz;i++){
		bm[i]=bv(i);
	}

	//3.4. Solve the Equaltions
	printf("SolveEqualtions");
	Eigen::VectorXd x = chol.solve(bv);
	double *result = new double[dimx*dimy*dimz];
	for(int i=0;i<dimx*dimy*dimz;i++){
		result[i]=x(i);
	}

	//4. Compare the Result
	//4.1. Get the Original Result
	double* originalResult = new double[dimx*dimy*dimz];
	double minResult=1000000, maxResult=-1000000;
	for(int i=0;i<dimz;i++)
	for(int j=0;j<dimy;j++)
	for(int k=0;k<dimx;k++){
		originalResult[i*dimx*dimy+j*dimx+k]=plusof(i,j,k,dimz,dimy,dimx);
		minResult = (minResult < originalResult[i*dimx*dimy+j*dimx+k])? minResult:originalResult[i*dimx*dimy+j*dimx+k];
		maxResult = (maxResult > originalResult[i*dimx*dimy+j*dimx+k])? maxResult:originalResult[i*dimx*dimy+j*dimx+k];
	}

	//4.2. Calculate the Errors
	cout<<"The resolution of the matrix is "<<dimx<<" x "<<dimy<<" x "<<dimz<<endl<<
		"In original result, max = "<<maxResult<<", min = "<<minResult<<endl<<
		"There are "<<gridNum<<"grid in calculation"<<endl;
	double error = absoluteError(result, originalResult, dimx*dimy*dimz);
	cout<<"error = "<<error<<endl<<
		"average error = "<<error/gridNum<<endl<<
		"relative error = "<<error/gridNum/(maxResult-minResult)/(maxResult-minResult)<<endl;

	//5. Show the result with OpenGL
	MatGraph3d * mat = new MatGraph3d(result, dimx, dimy, dimz, 0.5, true);
	mat->SetMask(poissonBuilder.GetHoleMask());
	SetMatGraph3d(mat);
	InitGL(argc,argv);
	return 0;
}
