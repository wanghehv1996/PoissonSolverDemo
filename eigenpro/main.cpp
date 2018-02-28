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



double absoluteError(double *mat1, double *mat2, int size){
	double totalerror = 0;
	for(int i=0;i<size;i++){
		totalerror += (mat1[i]-mat2[i]) * (mat1[i]-mat2[i]);
	}
	return totalerror;
}

int main(int argc, char** argv){

	//solve Ax = b; 
	
	//1. Build the problem
	int dimx = 50;
	int dimy = 50;
	int gridNum = dimx*dimy;
	PoissonBuilder poissonBuilder(dimy, dimx, SQRTPI*2/(dimx-1));
	poissonBuilder.SetForce(sinpowff);
	poissonBuilder.SetGradientX(sinpowgx);
	poissonBuilder.SetGradientY(sinpowgy);
	poissonBuilder.SetBoundary(PB_BOUND_Dirichlet | 
		PB_BOUND_X0 | PB_BOUND_X1 |
		PB_BOUND_Y0 | PB_BOUND_Y1,
		sinpowbc);
	// poissonBuilder.SetBoundary(PB_BOUND_Neumann | 
	// 	PB_BOUND_X0 | PB_BOUND_X1,
	// 	sinpowgx);
	// poissonBuilder.SetBoundary(PB_BOUND_Neumann | 
	// 	PB_BOUND_Y0 | PB_BOUND_Y1,
	// 	sinpowgy);

	//2. Set Some Holes
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -SQRTPI,SQRTPI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -SQRTPI,SQRTPI);
	for(int I=0;I<5;I++)
		for(int J=0;J<5;J++)
			for(int i=2;i<8;i++)
				for(int j=2;j<8;j++){
					double y = arry(I*10+i);
					double x = arrx(J*10+j);
					poissonBuilder.SetHole(I*10+i,J*10+j,PB_BOUND_Neumann,sin(y*y+x*x));
					//if set the hole, gridNum --
					gridNum--;
				}

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

	Eigen::MatrixXd bm(dimy,dimx);
	for(int i=0;i<dimy;i++)
		for(int j=0;j<dimx;j++){
			bm(i,j) = bv(i*dimx+j);
		}
	//cout<<endl<< "b = \n"<<bm<<endl;

	//3.4. Solve the Equaltions
	Eigen::VectorXd x = chol.solve(bv);
	double *result = new double[dimx*dimy];
	for(int i=0;i<dimx*dimy;i++){
		result[i]=x(i);
	}


	//4. Compare the Result
	//4.1. Get the Original Result
	double* originalResult = new double[dimx*dimy];
	double minResult=1000000, maxResult=-1000000;
	for(int i=0;i<dimy;i++){
		for(int j=0;j<dimx;j++){
			originalResult[i*dimx+j]=sinpowof(i,j,dimy,dimx);
			minResult = (minResult < originalResult[i*dimx+j])? minResult:originalResult[i*dimx+j];
			maxResult = (maxResult > originalResult[i*dimx+j])? maxResult:originalResult[i*dimx+j];
		}
	}

	//4.2. Calculate the Errors
	cout<<"The resolution of the matrix is "<<dimx<<" x "<<dimy<<endl<<
		"In original result, max = "<<maxResult<<", min = "<<minResult<<endl<<
		"There are "<<gridNum<<"grid in calculation"<<endl;
	double error = absoluteError(result, originalResult, dimx*dimy);
	cout<<"error = "<<error<<endl<<
		"average error = "<<error/gridNum<<endl<<
		"relative error = "<<error/gridNum/(maxResult-minResult)/(maxResult-minResult)<<endl;

	//5. Show the result with OpenGL
	MatGraph2d * mat = new MatGraph2d(result,dimx,dimy, 0.2,1,true);
	mat->SetMask(poissonBuilder.GetHoleMask());
	SetMatGraph2d(mat);
	InitGL(argc,argv);

	return 0;
}
