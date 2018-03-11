#include "PoissonBuilder.h"
#include <iostream>
using namespace std;

/******************************************************
 * Function:GetDebugBoundaryTypeStr(int boundaryinfo)
 *
 * return the debug string of boundary type
 ******************************************************/
char* GetDebugBoundaryTypeStr(int boundaryinfo){
	if((boundaryinfo & (PB_BOUND_Dirichlet | PB_BOUND_Neumann))==3 ||
		(boundaryinfo & (PB_BOUND_Dirichlet | PB_BOUND_Neumann))==0)
		return NULL;
	if(boundaryinfo & PB_BOUND_Dirichlet)
		return "Dirichlet";
	if(boundaryinfo & PB_BOUND_Neumann)
		return "Neumann";
	return NULL;
}

/******************************************************
 * Function:GetBoundaryPara(int boundarytype)
 *
 * If Dirichlet boundary condition. 
 *     ub - 2u0 + u1 = f(u0) && ub = sth 
 *     =>-2u0 + u1 + 0*u1 = f(u0) - sth
 *     so return 0
 *
 * If Neumann boundary condition.
 *     ub - 2u0 + u1 = f(u0) && u0' = sth ( so ub = u1 - sth * 2 * dx)
 *     =>-2u0 + u1 + 1*u1 = f(u0) + sth * 2 * dx
 *     so return 1
 ******************************************************/
int GetBoundaryPara(int boundarytype){
	if(boundarytype & PB_BOUND_Dirichlet)
		return 0;
	else if(boundarytype & PB_BOUND_Neumann)
		return 1;
}


int PoissonBuilder::SetForce(double (*forceFunc)(int, int, int, int)){
	_forceFunc = forceFunc;
}

int PoissonBuilder::SetGradientX(double (*gradientXFunc)(double, double, int, int)){
	_gradientXFunc = gradientXFunc;
}

int PoissonBuilder::SetGradientY(double (*gradientYFunc)(double, double, int, int)){
	_gradientYFunc = gradientYFunc;
}

/******************************************************
 * Function:SetBoundary(int boundarytype, double (*boundaryFunc)(int,int))
 * 
 * Class: PoissonBuilder
 * Set boundary type & function.
 * 
 ******************************************************/
int PoissonBuilder::SetBoundary(int boundarytype, double (*boundaryFunc)(double,double,int,int)){
	//x=0 boundary
	//store in _boundary[0]
	if(boundarytype & PB_BOUND_X0){
		_boundaryType[0] = boundarytype & 0x3;
		_boundaryFunc[0] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where x = 0 \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//x=dimx-1 boundary
	//store in _boundary[1]
	if(boundarytype & PB_BOUND_X1){
		_boundaryType[1] = boundarytype & 0x3;
		_boundaryFunc[1] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where x = colnum \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//y=0 boundary
	//store in _boundary[2]
	if(boundarytype & PB_BOUND_Y0){
		_boundaryType[2] = boundarytype & 0x3;
		_boundaryFunc[2] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where y = 0 \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//y=dimy-1 boundary
	//store in _boundary[3]
	if(boundarytype & PB_BOUND_Y1){
		_boundaryType[3] = boundarytype & 0x3;
		_boundaryFunc[3] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where y = rownum \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
}


int PoissonBuilder::SetHole(int y, int x, int type, double v){
	if(y<0 || y>=_dimy || x<0 || x>=_dimx)
		return -1;
	_holeMaskMap[y*_dimx+x]=type;
	_holeBoundaryMap[y*_dimx+x]=v;
}

int* PoissonBuilder::GetHoleMask(){
	return _holeMaskMap;
}

/******************************************************
 * Function:GetLaplacianMatrix()
 * 
 * Class:PoissonBuilder
 * 
 * build the laplacian matrix (A in Ax=b)
 * be careful about the boundary condition:
 *
 * If Dirichlet boundary condition. 
 *     ub - 2u0 + u1 = f(u0) && ub = sth 
 *     =>-2u0 + u1 + 0*u1 = f(u0) - sth
 *
 * If Neumann boundary condition.
 *     ub - 2u0 + u1 = f(u0) && (u0-ub)/dx = sth ( so ub = u0 - sth * dx)
 *     =>-2u0+u0 + u1 = f(u0) + sth * dx
 ******************************************************/
Eigen::SparseMatrix<double> PoissonBuilder::GetLaplacianMatrix(){

	int m = _dimx * _dimy;
	Eigen::SparseMatrix<double> A(m,m);
	std::vector<Eigen::Triplet<double> > coefficients;


	//build A
	for(int i=0;i<_dimy;i++)
		for(int j=0;j<_dimx;j++)
		{
			//skip if it's hole
			if(_holeMaskMap[i*_dimx+j]!=0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j, 1));
				continue;
			}

			int center = -4;

			//para of A(i+1,j) & A(i-1,j)
			//if y=0 boundary	
			if(i==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i+1)*_dimx+j, 1));
				center += (_boundaryType[2]==PB_BOUND_Neumann);
			}
			//if y=1 boundary
			else if(i==_dimy-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i-1)*_dimx+j, 1));
				center += (_boundaryType[3]==PB_BOUND_Neumann);
			}
			//no boundary in y-axis direction
			else{
				if(_holeMaskMap[(i+1)*_dimx+j]==0)//Normal
					coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i+1)*_dimx+j, 1));
				else if(_holeMaskMap[(i+1)*_dimx+j]==PB_BOUND_Neumann)//Neumann
					center += 1;
				
				if(_holeMaskMap[(i-1)*_dimx+j]==0)//Normal
					coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i-1)*_dimx+j, 1));
				else if(_holeMaskMap[(i-1)*_dimx+j]==PB_BOUND_Neumann)//Neumann
					center += 1;
				//printf("(%d,%d) ^%d v%d\t",i,j,iminus, iplus);
			}
			
			//para of A(i,j+1) & A(i,j-1)
			//if x=0 boundary
			if(j==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j+1, 1 ));
				center += (_boundaryType[0]==PB_BOUND_Neumann);
			}
			//if x=dimx-1 boundary
			else if(j==_dimx-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j-1, 1 ));
				center += (_boundaryType[1]==PB_BOUND_Neumann);
			}
			//no boundary in x-axis direction
			else{
				if(_holeMaskMap[i*_dimx+j+1]==0)//Normal
					coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j+1, 1 ));
				else if(_holeMaskMap[i*_dimx+j+1]==PB_BOUND_Neumann)//Neumann
					center += 1;
				
				if(_holeMaskMap[i*_dimx+j-1]==0)//Normal
					coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j-1, 1 ));
				else if(_holeMaskMap[i*_dimx+j-1]==PB_BOUND_Neumann)//Neumann
					center += 1;		

				
				//printf("(%d,%d) <%d >%d\t",i,j, jminus, jplus);
			}
			
			//para of A(i,j)
			coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j, center));
		}
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	return A;
}

/******************************************************
 * Function:GetForceVector()
 * 
 * Class:PoissonBuilder
 * 
 * build the force vector (b in Ax=b)
 * be careful about the boundary condition:
 *
 * If Dirichlet boundary condition. 
 *     ub - 2u0 + u1 = f(u0) && ub = sth 
 *     =>-2u0 + u1 + 0*u1 = f(u0) - sth
 *     so use force - sth
 *
 * If Neumann boundary condition.
 *     ub - 2u0 + u1 = f(u0) && u0' = sth ( so ub = u1 - sth * 2 * dx)
 *     =>-2u0 + u1 + 1*u1 = f(u0) + sth * 2 * dx
 *     so use force + sth * 2 * dx
 ******************************************************/
Eigen::VectorXd PoissonBuilder::GetForceVector(){
	Eigen::VectorXd b(_dimx*_dimy);

	//fill vector with force function
	for(int i=0;i<_dimy;i++){
		for(int j=0;j<_dimx;j++)
		{	
			if(_holeMaskMap[i*_dimx+j] == 0)
				b(i*_dimx+j)=_forceFunc(i,j,_dimy,_dimx)*_dx*_dx;
		}
	}

	int i,j;

	//use different boundary condition
	j=0;
	//x=0 boundary
	for(i = 0; i < _dimy; i++){
		if(_holeMaskMap[i*_dimx+j]==0)
		if(_boundaryType[0] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (-1,i)
			b(i*_dimx+j)-=_boundaryFunc[0](i,-1, _dimy, _dimx);
		else if(_boundaryType[0] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (-0.5,i)
			b(i*_dimx+j)+=_boundaryFunc[0](i, -0.5, _dimy, _dimx) * _dx;
	}

	//x=dimx-1 boundary
	j=_dimx-1;
	for(i = 0; i < _dimy; i++){
		if(_holeMaskMap[i*_dimx+j]==0)
		if(_boundaryType[1] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (dimx,i)
			b(i*_dimx+j)-=_boundaryFunc[1](i, _dimx, _dimy, _dimx);
		else if(_boundaryType[1] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (dimx-0.5,i)
			b(i*_dimx+j)-=_boundaryFunc[1](i, _dimx-0.5, _dimy, _dimx) * _dx;
	}

	//y=0 boundary
	i=0;
	for(j = 0; j < _dimx; j++){
		if(_holeMaskMap[i*_dimx+j]==0)
		if(_boundaryType[2] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (j,-1)
			b(i*_dimx+j)-=_boundaryFunc[2](-1, j, _dimy, _dimx);
		else if(_boundaryType[2] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (j,-0.5)
			b(i*_dimx+j)+=_boundaryFunc[2](-0.5, j, _dimy, _dimx) * _dx;
	}

	//y=dimy-1 boundary
	i=_dimy-1;
	for(j = 0; j < _dimx; j++){
		if(_holeMaskMap[i*_dimx+j]==0)
		if(_boundaryType[3] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (j,dimy)
			b(i*_dimx+j)-=_boundaryFunc[3](_dimy, j, _dimy, _dimx);
		else if(_boundaryType[3] & PB_BOUND_Neumann)
			//Neumann, boundary is the value on (j,dimy-1)
			b(i*_dimx+j)-=_boundaryFunc[3](_dimy-0.5, j, _dimy, _dimx) * _dx ;
	}

	for(i = 0; i < _dimy; i++)
		for(j = 0; j < _dimx; j++)
			if(_holeMaskMap[i*_dimx+j] == 0){//not hole
				//y-1 neighbour is hole
				if(i-1>0 && _holeMaskMap[(i-1)*_dimx+j])

					if(_holeMaskMap[(i-1)*_dimx+j]==1)
						b(i*_dimx+j)-=_holeBoundaryMap[(i-1)*_dimx+j];
					else if(_holeMaskMap[(i-1)*_dimx+j]==2)
						b(i*_dimx+j)+= _dx * _gradientYFunc(i-0.5,j,_dimy,_dimx);

				//y+1 neighbour is hole
				if(i+1<_dimy && _holeMaskMap[(i+1)*_dimx+j])

					if(_holeMaskMap[(i+1)*_dimx+j]==1)
						b(i*_dimx+j)-=_holeBoundaryMap[(i+1)*_dimx+j];
					else if(_holeMaskMap[(i+1)*_dimx+j]==2)
						b(i*_dimx+j)-= _dx * _gradientYFunc(i+0.5,j,_dimy,_dimx);

				//x-1 neighbour is hole
				if(j-1>0 && _holeMaskMap[i*_dimx+j-1])

					if(_holeMaskMap[i*_dimx+j-1]==1)
						b(i*_dimx+j)-=_holeBoundaryMap[i*_dimx+j-1];
					else if(_holeMaskMap[i*_dimx+j-1]==2)
						b(i*_dimx+j)+= _dx * _gradientXFunc(i,j-0.5,_dimy,_dimx);

				//x+1 neighbour is hole
				if(j+1<_dimx && _holeMaskMap[i*_dimx+j+1])

					if(_holeMaskMap[i*_dimx+j+1]==1)
						b(i*_dimx+j)-=_holeBoundaryMap[i*_dimx+j+1];
					else if(_holeMaskMap[i*_dimx+j+1]==2)
						b(i*_dimx+j)-= _dx * _gradientXFunc(i,j+0.5,_dimy,_dimx);
			}else{//is a hole
				b(i*_dimx+j) = sinmulof(i,j,_dimy,_dimx);//_holeBoundaryMap[i*_dimx+j];
				//cout<<"leak b = "<<_holeBoundaryMap[i*_dimx+j]<<endl;
			}

	return b;
}

double zerof(double yi, double xi, int dimy, int dimx){
	return 0;
}

double sinpowof(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -SQRTPI,SQRTPI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -SQRTPI,SQRTPI);
	double x = arrx(xi);
	double y = arry(yi);
	return sin(x*x+y*y);
}

//force function of 
//z = sin(x*x+y*y)
//laplacian(z) = 4*cos(x*x+y*y) - 4*sin(x*x+y*y)*(x*x+y*y)
double sinpowff(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -SQRTPI,SQRTPI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -SQRTPI,SQRTPI);
	double x = arrx(xi);
	double y = arry(yi);
	return 4*cos(x*x + y*y) - sin(x*x + y*y)*4*(x*x+y*y);

}

//boundary condition of
//z = sin(x*x+y*y)
double sinpowbc(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -SQRTPI+xi*2*SQRTPI/(dimx-1);
	y = -SQRTPI+yi*2*SQRTPI/(dimy-1);

	
	return sin(x*x+y*y);
}

//gradient condition on x-axis of 
//z = sin(x*x+y*y)
//dz/dx = 2 * x * cos(x*x+y*y)
double sinpowgx(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -SQRTPI+xi*2*SQRTPI/(dimx-1);
	y = -SQRTPI+yi*2*SQRTPI/(dimy-1);
	return x * 2 * cos(x*x+y*y);
}

//gradient condition on y-axis of 
//z = sin(x*x+y*y)
//dz/dy = 2 * y * cos(x*x+y*y)
double sinpowgy(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -SQRTPI+xi*2*SQRTPI/(dimx-1);
	y = -SQRTPI+yi*2*SQRTPI/(dimy-1);
	return y * 2 * cos(x*x+y*y);
}

//---------------------------------------------
//z = sin(x) * sin(y)
double sinmulof(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -PI,PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -PI,PI);
	double x = arrx(xi);
	double y = arry(yi);
	return sin(x)*sin(y);
}

//force function of 
//z = sin(x) * sin(y)
//laplacian(z) = -2*sin(x)*sin(y)
double sinmulff(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -PI,PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -PI,PI);
	double x = arrx(xi);
	double y = arry(yi);
	return -2*sin(x)*sin(y);

}

//boundary condition of
//z = sin(x) * sin(y)
double sinmulbc(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return sin(x)*sin(y);
}

//gradient condition on x-axis of 
//z = sin(x) * sin(y)
//dz/dx =  cos(x) * sin(y)
double sinmulgx(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return cos(x) * sin(y);
}

//gradient condition on y-axis of 
//z = sin(x) * sin(y)
//dz/dy = sin(x) * cos(y)
double sinmulgy(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return sin(x) * cos(y);
}
//-------------------------------------------------

//-------------------------------------------------
//z = sin(2x) * sin(2y)
double sinmul2of(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -PI,PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -PI,PI);
	double x = arrx(xi);
	double y = arry(yi);
	return sin(2*x)*sin(2*y);
}

//force function of 
//z = sin(2x) * sin(2y)
//laplacian(z) = -8*sin(x)*sin(y)
double sinmul2ff(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -PI,PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -PI,PI);
	double x = arrx(xi);
	double y = arry(yi);
	return -8*sin(2*x)*sin(2*y);

}

//boundary condition of
//z = sin(2x) * sin(2y)
double sinmul2bc(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return sin(2*x)*sin(2*y);
}

//gradient condition on x-axis of 
//z = sin(2x) * sin(2y)
//dz/dx = 2*cos(2x) * sin(2y)
double sinmul2gx(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return 2 * cos(2*x) * sin(2*y);
}

//gradient condition on y-axis of 
//z = sin(x) * sin(y)
//dz/dy = 2*sin(2x) * cos(2y)
double sinmul2gy(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return 2 * sin(2*x) * cos(2*y);
}
//-------------------------------------------------

//-------------------------------------------------
//z = sin(4x) * sin(4y)
double sinmul4of(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -PI,PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -PI,PI);
	double x = arrx(xi);
	double y = arry(yi);
	return sin(4*x)*sin(4*y);
}

//force function of 
//z = sin(4x) * sin(4y)
//laplacian(z) = -8*sin(4x)*sin(4y)
double sinmul4ff(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -PI,PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -PI,PI);
	double x = arrx(xi);
	double y = arry(yi);
	return -32*sin(4*x)*sin(4*y);

}

//boundary condition of
//z = sin(4x) * sin(4y)
double sinmul4bc(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return sin(4*x)*sin(4*y);
}

//gradient condition on x-axis of 
//z = sin(4x) * sin(4y)
//dz/dx = 4*cos(4x) * sin(4y)
double sinmul4gx(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return 4 * cos(4*x) * sin(4*y);
}

//gradient condition on y-axis of 
//z = sin(4x) * sin(4y)
//dz/dy = 4*sin(4x) * cos(4y)
double sinmul4gy(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return 4 * sin(4*x) * cos(4*y);
}
//-------------------------------------------------

//-------------------------------------------------
//z = sin(8*x) * sin(8*y)
double sinmul8of(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -PI,PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -PI,PI);
	double x = arrx(xi);
	double y = arry(yi);
	return sin(8*x)*sin(8*y);
}

//force function of 
//z = sin(8x) * sin(8y)
//laplacian(z) = -8*sin(8x)*sin(8y)
double sinmul8ff(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -PI,PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -PI,PI);
	double x = arrx(xi);
	double y = arry(yi);
	return -128*sin(8*x)*sin(8*y);

}

//boundary condition of
//z = sin(8x) * sin(8y)
double sinmul8bc(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return sin(8*x)*sin(8*y);
}

//gradient condition on x-axis of 
//z = sin(4x) * sin(4y)
//dz/dx = 4*cos(4x) * sin(4y)
double sinmul8gx(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return 8 * cos(8*x) * sin(8*y);
}

//gradient condition on y-axis of 
//z = sin(8x) * sin(8y)
//dz/dy = 8*sin(8x) * cos(8y)
double sinmul8gy(double yi, double xi, int dimy, int dimx){
	double x,y;	
	x = -PI+xi*2*PI/(dimx-1);
	y = -PI+yi*2*PI/(dimy-1);
	return 8 * sin(8*x) * cos(8*y);
}
//-------------------------------------------------

//z = x+y
double plusff(int yi, int xi, int dimy, int dimx){
	return 0;
}


double plusbc(int yi, int xi, int dimy, int dimx){
	return yi+xi;
}

double plusgx(int yi, int xi, int dimy, int dimx){
	return 1;
}
double plusgy(int yi, int xi, int dimy, int dimx){
	return 1;
}