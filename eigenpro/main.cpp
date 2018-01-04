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

double sinforce(int r,int c, int i, int j)
{
	Eigen::ArrayXd boundaryr = Eigen::ArrayXd::LinSpaced(r, 0,M_PI).sin().pow(2);

	Eigen::ArrayXd boundaryc = Eigen::ArrayXd::LinSpaced(c, 0,M_PI).sin().pow(2);
	double result=0;
	if(i==0||i==r-1)
		result += boundaryc(j);
	if(j==0||j==c-1)
		result += boundaryr(i);
	return result;
}

double linearforce(int r,int c, int i, int j)
{
	Eigen::ArrayXd boundaryr = Eigen::ArrayXd::LinSpaced(r, 0,1);

	Eigen::ArrayXd boundaryc = Eigen::ArrayXd::LinSpaced(c, 0,1);
	double result = boundaryc(j) + boundaryr(i);
	return result;
}

double zeroforce(int r,int c, int i, int j)
{
	return 0;	
}

double sinpowboundary(int i, int total){
	ArrayXd arr= Eigen::ArrayXd::LinSpaced(total, 0,M_PI).sin().pow(2);
	return arr(i);
}

double cospowboundary(int i, int total){
	ArrayXd arr= Eigen::ArrayXd::LinSpaced(total, 0,M_PI).cos().pow(2);
	return arr(i);
}

double parabolapowboundary(int i, int total){
	ArrayXd arr= Eigen::ArrayXd::LinSpaced(total, -1,1).pow(2);
	return arr(i);
}

//Dirichlet boundary condition
VectorXd AddDirichletBC(int r, int c, VectorXd& b, double (*boundaryFunc)(int, int)){
	int i,j;

	i=0;
	for(j = 0; j < c; j++){
		b(i*r+j)-=boundaryFunc(j, c);
	}
	i=r-1;
	for(j = 0; j < c; j++){
		b(i*r+j)-=boundaryFunc(j, c);
	}

	j=0;
	for(i = 0; i < r; i++){
		b(i*r+j)-=boundaryFunc(i, r);
	}
	j=c-1;
	for(i = 0; i < r; i++){
		b(i*r+j)-=boundaryFunc(i, r);
	}
	return b;
}

//Neumann boundary condition
VectorXd AddNeumannBC(int r, int c, VectorXd& b, double (*boundaryFunc)(int, int), double dx = 1){
	int i,j;

	i=0;
	for(j = 0; j < c; j++){
		b(i*r+j)+=boundaryFunc(j, c) * dx * 2;
	}
	i=r-1;
	for(j = 0; j < c; j++){
		b(i*r+j)+=boundaryFunc(j, c) * dx * 2;
	}

	j=0;
	for(i = 0; i < r; i++){
		b(i*r+j)+=boundaryFunc(i, r) * dx * 2;
	}
	j=c-1;
	for(i = 0; i < r; i++){
		b(i*r+j)+=boundaryFunc(i, r) * dx * 2;
	}
	return b;
}

VectorXd GetForce(int r, int c, double (*forceFunc)(int,int,int,int)){
	VectorXd force(r*c);
	for(int i=0;i<r;i++){
		for(int j=0;j<c;j++)
		{	
			force(i*r+j)=forceFunc(r,c,i,j);
		}
	}
	return force;
}

SparseMatrix<double> buildLaplacianMatrix(int r, int c, int boundStrategy = 0)
{
	int m = r*c;
	SparseMatrix<double> A(m,m);
	vector<Triplet<double> > coefficients;


	//build A
	for(int i=0;i<r;i++)
		for(int j=0;j<c;j++)
		{
			//if boundStrategy=0, it is Dirichlet boundary condition. 
			//    ub - 2u0 + u1 = f(u0) && ub = sth 
			//  =>-2u0 + u1 = f(u0) - sth

			//if boundStrategy=1, it is Neumann boundary condition.
			//    ub - 2u0 + u1 = f(u0) && u0' = sth ( so ub = u1 - sth * 2 * dx)
			//  =>-2u0 + 2u1 = f(u0) + sth * 2 * dx
			
			if(i==0){
				coefficients.push_back(Triplet<double>( i*c+j, (i+1)*c+j, 1 + boundStrategy));
			}
			else if(i==r-1){
				coefficients.push_back(Triplet<double>( i*c+j, (i-1)*c+j, 1 + boundStrategy));
			}
			else{
				coefficients.push_back(Triplet<double>( i*c+j, (i+1)*c+j, 1));
				coefficients.push_back(Triplet<double>( i*c+j, (i-1)*c+j, 1));
			}
				
			if(j==0){
				coefficients.push_back(Triplet<double>( i*c+j, i*c+j+1, 1 + boundStrategy));
			}else if(j==c-1){
				coefficients.push_back(Triplet<double>( i*c+j, i*c+j-1, 1 + boundStrategy));
			}else{

				coefficients.push_back(Triplet<double>( i*c+j, i*c+j+1, 1));
				coefficients.push_back(Triplet<double>( i*c+j, i*c+j-1, 1));
			}
			
			coefficients.push_back(Triplet<double>( i*c+j, i*c+j, -4));//center));
		}
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	return A;
}

//laplacian(field) * x = b
MatrixXd solveLalacianProblem(int r, int c, double (*forceFunc)(int,int,int,int), int boundStrategy = 0){

	SparseMatrix<double> A = buildLaplacianMatrix(r, c, boundStrategy);
	//show A
	cout<<endl<<"A = \n"<<A<<endl;

	//build b
	VectorXd bv = GetForce(r,c,forceFunc);

	//use different boundary condition
	if(boundStrategy == 0)
		bv = AddDirichletBC(r,c,bv,cospowboundary);//parabolapowboundary, cospowboundary
	else if(boundStrategy == 1)
		bv = AddNeumannBC(r,c,bv,sinpowboundary);

	//show b
	Eigen::MatrixXd bm(r,c);
	for(int i=0;i<r;i++)
		for(int j=0;j<c;j++){
			bm(i,j) = bv(i*r+j);
		}
	cout<<endl<< "b = \n"<<bm<<endl;

	//solve Ax=b
	//you can use different solver here
	Eigen::SimplicialCholesky<SparseMatrix<double> > chol(A);// performs a Cholesky factorization of A
	Eigen::VectorXd x = chol.solve(bv);// use the factorization to solve for the given right hand side

	//return x
	MatrixXd xm(r,c);
	for(int i=0;i<r;i++)
		for(int j=0;j<c;j++){
			xm(i,j) = x(i*c+j);
		}
	return xm;
}

int main(int argc, char** argv){
	int r=15,c=15;

	//solve Ax = b; 
	//arg:row num, col num, forcefunction(b)
	MatrixXd result = solveLalacianProblem(r,c,zeroforce,0);
	cout<<endl<<"result = \n"<<result<<endl;

	//arg: matrix, horizontal scale, vertical scale, mesh or line?
	MatGraph * mat = new MatGraph(result, 0.2,0.5,true);
	SetMatGraph(mat);
	InitGL(argc,argv);
	return 0;
}
