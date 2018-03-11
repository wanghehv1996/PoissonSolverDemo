#include "eigen3/Eigen/Dense"
#include <eigen3/Eigen/Eigen> 
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <time.h>
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

double levelError(double *mat1, double *mat2, int dimy, int dimx){
	double totalerror = 0;
	int dimBlock = 10;
	int Imax = dimy/dimBlock;
	int Jmax = dimx/dimBlock;

	for(int level = 0;level<dimBlock/2;level++){
		totalerror = 0;
		for(int I=0;I<Imax;I++)
		for(int J=0;J<Jmax;J++){

			for(int j = level;j<dimBlock-level;j++){
				int index = (I*dimBlock+level)*dimx + J*dimBlock+j;
				//cout<<'('<<index/dimx<<' '<<index%dimx<<')';
				totalerror += (mat1[index]-mat2[index]) * (mat1[index]-mat2[index]);
			}
			for(int i = level+1;i<dimBlock-level-1;i++){
				int index = (I*dimBlock+i)*dimx + J*dimBlock+level;
				//cout<<'('<<index/dimx<<' '<<index%dimx<<')';
				totalerror += (mat1[index]-mat2[index]) * (mat1[index]-mat2[index]);
				index = (I*dimBlock+i)*dimx + J*dimBlock+dimBlock-level-1;
				//cout<<'('<<index/dimx<<' '<<index%dimx<<')';
				totalerror += (mat1[index]-mat2[index]) * (mat1[index]-mat2[index]);
			}
			for(int j = level;j<dimBlock-level;j++){
				int index = (I*dimBlock+dimBlock-level-1)*dimx + J*dimBlock+j;
				//cout<<'('<<index/dimx<<' '<<index%dimx<<')';
				totalerror += (mat1[index]-mat2[index]) * (mat1[index]-mat2[index]);
			}
		}
		cout//<<"Level "<<level<<endl<<
			//"error = "<<totalerror<<endl<<
			//"error_avg = "<<totalerror/(Imax*Jmax*(dimBlock-level-1)*4)<<endl<<
			//"sqrt error_avg = "	
			<< sqrt(totalerror/(Imax*Jmax*(dimBlock-level-1)*4))<<//endl<<
			//"gridNum = "<<(Imax*Jmax*(dimBlock-2*level-1)*4)<<endl<<
			endl;
	}

	return 0;
}

int main(int argc, char** argv){

	//solve Ax = b; 
	
	//1. Build the problem
	int dimx = 50;
	int dimy = 50;
	int gridNum = dimx*dimy;
	// PoissonBuilder poissonBuilder(dimy, dimx, SQRTPI*2/(dimx-1));
	// poissonBuilder.SetForce(sinpowff);
	// poissonBuilder.SetGradientX(sinpowgx);
	// poissonBuilder.SetGradientY(sinpowgy);
	// poissonBuilder.SetBoundary(PB_BOUND_Dirichlet | 
	// 	PB_BOUND_X0 | PB_BOUND_X1 |
	// 	PB_BOUND_Y0 | PB_BOUND_Y1,
	// 	sinpowbc);
	// // poissonBuilder.SetBoundary(PB_BOUND_Neumann | 
	// // 	PB_BOUND_X0 | PB_BOUND_X1,
	// // 	sinpowgx);
	// // poissonBuilder.SetBoundary(PB_BOUND_Neumann | 
	// // 	PB_BOUND_Y0 | PB_BOUND_Y1,
	// // 	sinpowgy);

	PoissonBuilder poissonBuilder(dimy, dimx, PI*2/(dimx-1));
	poissonBuilder.SetForce(sinmulff);
	poissonBuilder.SetGradientX(zerof);//sinmulgx);
	poissonBuilder.SetGradientY(zerof);//sinmulgy);//zerof);//
	poissonBuilder.SetBoundary(PB_BOUND_Dirichlet | 
		PB_BOUND_X0 | PB_BOUND_X1 |
		PB_BOUND_Y0 | PB_BOUND_Y1,
		sinmulbc);

	//2. Set Some Holes
	//Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -SQRTPI,SQRTPI);
	//Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -SQRTPI,SQRTPI);
	for(int I=0;I<5;I++)
		for(int J=0;J<5;J++)
			for(int i=2;i<8;i++)
				for(int j=2;j<8;j++){
					//double y = arry(I*10+i);
					//double x = arrx(J*10+j);
					int ii = I*10+i;
					int jj = J*10+j;
					poissonBuilder.SetHole(ii,jj,PB_BOUND_Dirichlet,0);//sinmulof(ii,jj,dimy,dimx));
					//if set the hole, gridNum --
					gridNum--;
				}

	//3. Solve the Problem Ax=b
	//3.1. Get Laplacian Matrix A
	printf("GetLaplacianMatrix\n");
	SparseMatrix<double> A = poissonBuilder.GetLaplacianMatrix();
	//cout<<endl<<"A1 = \n"<<A<<endl;

	//3.2. Try Different Solver
	Eigen::ConjugateGradient<SparseMatrix<double> > chol(A);
	
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
	long time_start = clock();
	Eigen::VectorXd x = chol.solve(bv);
	long time_end = clock();
	printf("Solve Equation\n");
	cout<<"Use Time:"<<(time_end - time_start) * 1.0 / CLOCKS_PER_SEC<<"s"<<endl;
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
			originalResult[i*dimx+j]=sinmulof(i,j,dimy,dimx);
			minResult = (minResult < originalResult[i*dimx+j])? minResult:originalResult[i*dimx+j];
			maxResult = (maxResult > originalResult[i*dimx+j])? maxResult:originalResult[i*dimx+j];
		}
	}

	double minResultF=1000000, maxResultF=-1000000;
	for(int i=0;i<dimy;i++){
		for(int j=0;j<dimx;j++){
			minResultF = (minResultF < result[i*dimx+j])? minResultF:result[i*dimx+j];
			maxResultF = (maxResultF > result[i*dimx+j])? maxResultF:result[i*dimx+j];
		}
	}
	cout<<"f_max = "<<maxResultF<<endl<<
		"f_min = "<<minResultF<<endl;

	//4.2. Calculate the Errors
	cout<<"Resolution "<<dimx<<" x "<<dimy<<endl<<
		"n_grid = "<<gridNum<<endl<<
		"phi_max = "<<maxResult<<", phi_min = "<<minResult<<endl;
	double error = absoluteError(result, originalResult, dimx*dimy);
	cout<<"error_total = "<<error<<endl<<
		"error_avg = "<<error/gridNum<<endl<<
		"error_rel = "<<error/gridNum/(maxResult-minResult)/(maxResult-minResult)<<endl;
	//levelError(result, originalResult, dimy, dimx);

	//4.3. Show Error Map
	double *errorResult = new double[dimx*dimy];
	for(int i=0;i<dimy;i++){
		for(int j=0;j<dimx;j++){
			int index = i*dimx + j;
			errorResult[index]=abs(result[index]-originalResult[index]);
		}
	}

	//5. Show the result with OpenGL
	MatGraph2d * mat = new MatGraph2d(result,dimx,dimy, 0.2,1,true);
	//mat->SetMask(poissonBuilder.GetHoleMask());
	SetMatGraph2d(mat);
	InitGL(argc,argv);

	return 0;
}
