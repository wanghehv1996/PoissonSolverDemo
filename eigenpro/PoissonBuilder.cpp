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

int PoissonBuilder::SetGradientX(double (*gradientXFunc)(int, int, int, int)){
	_gradientXFunc = gradientXFunc;
}

int PoissonBuilder::SetGradientY(double (*gradientYFunc)(int, int, int, int)){
	_gradientYFunc = gradientYFunc;
}

/******************************************************
 * Function:SetBoundary(int boundarytype, double (*boundaryFunc)(int,int))
 * 
 * Class: PoissonBuilder
 * Set boundary type & function.
 * 
 ******************************************************/
int PoissonBuilder::SetBoundary(int boundarytype, double (*boundaryFunc)(int,int,int,int)){
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
 *     ub - 2u0 + u1 = f(u0) && u0' = sth ( so ub = u1 - sth * 2 * dx)
 *     =>-2u0 + u1 + 1*u1 = f(u0) + sth * 2 * dx
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


			//para of A(i+1,j) & A(i-1,j)
			//if y=0 boundary	
			if(i==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i+1)*_dimx+j, 1 + GetBoundaryPara(_boundaryType[2])));
			}
			//if y=1 boundary
			else if(i==_dimy-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i-1)*_dimx+j, 1 + GetBoundaryPara(_boundaryType[3])));
			}
			//no boundary in y-axis direction
			else{
				int iplus = 0;
				int iminus = 0;

				if(_holeMaskMap[(i+1)*_dimx+j]==0)//Normal
					iplus += 1;
				else if(_holeMaskMap[(i+1)*_dimx+j]==1)//Dirichlet
					;//iplus = 0;
				else if(_holeMaskMap[(i+1)*_dimx+j]==2)//Neumann
					iminus += 1;
				
				if(_holeMaskMap[(i-1)*_dimx+j]==0)//Normal
					iminus += 1;
				else if(_holeMaskMap[(i-1)*_dimx+j]==1)//Dirichlet
					;//iminus = 0;
				else if(_holeMaskMap[(i-1)*_dimx+j]==2)//Neumann
					iplus += 1;				

				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i+1)*_dimx+j, iplus));
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i-1)*_dimx+j, iminus));

				//printf("(%d,%d) ^%d v%d\t",i,j,iminus, iplus);
			}
			
			//para of A(i,j+1) & A(i,j-1)
			//if x=0 boundary
			if(j==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j+1, 1 + GetBoundaryPara(_boundaryType[0]) ));
			}
			//if x=dimx-1 boundary
			else if(j==_dimx-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j-1, 1 + GetBoundaryPara(_boundaryType[1]) ));
			}
			//no boundary in x-axis direction
			else{
				int jplus = 0;
				int jminus = 0;

				if(_holeMaskMap[i*_dimx+j+1]==0)//Normal
					jplus += 1;
				else if(_holeMaskMap[i*_dimx+j+1]==1)//Dirichlet
					;//jplus = 0;
				else if(_holeMaskMap[i*_dimx+j+1]==2)//Neumann
					jminus += 1;
				
				if(_holeMaskMap[i*_dimx+j-1]==0)//Normal
					jminus += 1;
				else if(_holeMaskMap[i*_dimx+j-1]==1)//Dirichlet
					;//jminus = 0;
				else if(_holeMaskMap[i*_dimx+j-1]==2)//Neumann
					jplus += 1;				

				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j+1, jplus));
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j-1, jminus));
				//printf("(%d,%d) <%d >%d\t",i,j, jminus, jplus);
			}
			
			//para of A(i,j)
			coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j, -4));//center));
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
			//Neumann, boundary is the gradient on (0,i)
			b(i*_dimx+j)+=_boundaryFunc[0](i, 0, _dimy, _dimx) * _dx * 2;
	}

	//x=dimx-1 boundary
	j=_dimx-1;
	for(i = 0; i < _dimy; i++){
		if(_holeMaskMap[i*_dimx+j]==0)
		if(_boundaryType[1] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (dimx,i)
			b(i*_dimx+j)-=_boundaryFunc[1](i, _dimx, _dimy, _dimx);
		else if(_boundaryType[1] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (dimx-1,i)
			b(i*_dimx+j)-=_boundaryFunc[1](i, _dimx-1, _dimy, _dimx) * _dx * 2;
	}

	//y=0 boundary
	i=0;
	for(j = 0; j < _dimx; j++){
		if(_holeMaskMap[i*_dimx+j]==0)
		if(_boundaryType[2] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (j,-1)
			b(i*_dimx+j)-=_boundaryFunc[2](-1, j, _dimy, _dimx);
		else if(_boundaryType[2] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (j,0)
			b(i*_dimx+j)+=_boundaryFunc[2](0, j, _dimy, _dimx) * _dx * 2;
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
			b(i*_dimx+j)-=_boundaryFunc[3](_dimy-1, j, _dimy, _dimx) * _dx * 2;
	}

	for(i = 0; i < _dimy; i++)
		for(j = 0; j < _dimx; j++)
			if(_holeMaskMap[i*_dimx+j] == 0){//not hole
				//y-1 neighbour is hole
				if(i-1>0 && _holeMaskMap[(i-1)*_dimx+j])

					if(_holeMaskMap[(i-1)*_dimx+j]==1)
						b(i*_dimx+j)-=_holeBoundaryMap[(i-1)*_dimx+j];
					else if(_holeMaskMap[(i-1)*_dimx+j]==2)
						b(i*_dimx+j)+= _dx * 2 * _gradientYFunc(i,j,_dimy,_dimx);

				//y+1 neighbour is hole
				if(i+1<_dimy && _holeMaskMap[(i+1)*_dimx+j])

					if(_holeMaskMap[(i+1)*_dimx+j]==1)
						b(i*_dimx+j)-=_holeBoundaryMap[(i+1)*_dimx+j];
					else if(_holeMaskMap[(i+1)*_dimx+j]==2)
						b(i*_dimx+j)-= _dx * 2 * _gradientYFunc(i,j,_dimy,_dimx);

				//x-1 neighbour is hole
				if(j-1>0 && _holeMaskMap[i*_dimx+j-1])

					if(_holeMaskMap[i*_dimx+j-1]==1)
						b(i*_dimx+j)-=_holeBoundaryMap[i*_dimx+j-1];
					else if(_holeMaskMap[i*_dimx+j-1]==2)
						b(i*_dimx+j)+= _dx * 2 * _gradientXFunc(i,j,_dimy,_dimx);

				//x+1 neighbour is hole
				if(j+1<_dimx && _holeMaskMap[i*_dimx+j+1])

					if(_holeMaskMap[i*_dimx+j+1]==1)
						b(i*_dimx+j)-=_holeBoundaryMap[i*_dimx+j+1];
					else if(_holeMaskMap[i*_dimx+j+1]==2)
						b(i*_dimx+j)-= _dx * 2 * _gradientXFunc(i,j,_dimy,_dimx);
			}else{//is a hole
				b(i*_dimx+j) = _holeBoundaryMap[i*_dimx+j];
				//cout<<"leak b = "<<_holeBoundaryMap[i*_dimx+j]<<endl;
			}

	return b;
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
double sinpowbc(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -SQRTPI,SQRTPI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -SQRTPI,SQRTPI);
	
	double x,y;
	if(xi>=0&&xi<dimx)
		x = arrx(xi);
	else
		x = SQRTPI*(dimx+1)/(dimx-1);
	
	if(yi>=0&&yi<dimy)
		y = arry(yi);
	else
		y = SQRTPI*(dimy+1)/(dimy-1);
	
	return sin(x*x+y*y);
	//printf("%f\t",arr(i));
	//return sin(arr(i)+(SQRTPI*(total+1)/(total-1))*(SQRTPI*(total+1)/(total-1)));
}

//gradient condition on x-axis of 
//z = sin(x*x+y*y)
//dz/dx = 2 * x * cos(x*x+y*y)
double sinpowgx(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -SQRTPI,SQRTPI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -SQRTPI,SQRTPI);
	double x = arrx(xi);
	double y = arry(yi);
	return x * 2 * cos(x*x+y*y);
}

//gradient condition on y-axis of 
//z = sin(x*x+y*y)
//dz/dy = 2 * y * cos(x*x+y*y)
double sinpowgy(int yi, int xi, int dimy, int dimx){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, -SQRTPI,SQRTPI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, -SQRTPI,SQRTPI);
	double x = arrx(xi);
	double y = arry(yi);
	return y * 2 * cos(x*x+y*y);
}

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