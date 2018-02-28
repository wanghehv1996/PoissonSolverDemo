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


int PoissonBuilder::SetForce(double (*forceFunc)(int, int, int, int, int, int)){
	_forceFunc = forceFunc;
}

int PoissonBuilder::SetGradientX(double (*gradientXFunc)(int, int, int, int, int, int)){
	_gradientXFunc = gradientXFunc;
}

int PoissonBuilder::SetGradientY(double (*gradientYFunc)(int, int, int, int, int, int)){
	_gradientYFunc = gradientYFunc;
}

int PoissonBuilder::SetGradientZ(double (*gradientZFunc)(int, int, int, int, int, int)){
	_gradientZFunc = gradientZFunc;
}

/******************************************************
 * Function:SetBoundary(int boundarytype, double (*boundaryFunc)(int,int))
 * 
 * Class: PoissonBuilder
 * Set boundary type & function.
 * 
 ******************************************************/
int PoissonBuilder::SetBoundary(int boundarytype, double (*boundaryFunc)(int,int,int,int,int,int)){
	//left boundary
	//store in _boundary[0]
	if(boundarytype & PB_BOUND_X0){
		_boundaryType[0] = boundarytype & 0x3;
		_boundaryFunc[0] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where x = 0 \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//right boundary
	//store in _boundary[1]
	if(boundarytype & PB_BOUND_X1){
		_boundaryType[1] = boundarytype & 0x3;
		_boundaryFunc[1] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where x = resx \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//up boundary
	//store in _boundary[2]
	if(boundarytype & PB_BOUND_Y0){
		_boundaryType[2] = boundarytype & 0x3;
		_boundaryFunc[2] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where y = 0 \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//down boundary
	//store in _boundary[3]
	if(boundarytype & PB_BOUND_Y1){
		_boundaryType[3] = boundarytype & 0x3;
		_boundaryFunc[3] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where y = resy \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//front boundary
	//store in _boundary[4]
	if(boundarytype & PB_BOUND_Z0){
		_boundaryType[4] = boundarytype & 0x3;
		_boundaryFunc[4] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where z = 0 \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//back boundary
	//store in _boundary[5]
	if(boundarytype & PB_BOUND_Z1){
		_boundaryType[5] = boundarytype & 0x3;
		_boundaryFunc[5] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where z = resz \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
}


int PoissonBuilder::SetHole(int z,int y, int x, int type, double v){
	if(z<0 || z>=_dimz || y<0 || y>=_dimy || x<0 || x>=_dimx)
		return -1;
	_holeMaskMap[z*_dimx*_dimy + y*_dimx + x]=type;
	_holeBoundaryMap[z*_dimx*_dimy + y*_dimx + x]=v;
	return 0;
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

	int m = _dimx * _dimy * _dimz;
	int _dimxy = _dimx * _dimy;
	Eigen::SparseMatrix<double> A(m,m);
	std::vector<Eigen::Triplet<double> > coefficients;

	//build A
	for(int i=0;i<_dimz;i++)
		for(int j=0;j<_dimy;j++)
			for(int k=0;k<_dimx;k++)
		{
			//skip if it's hole
			if(_holeMaskMap[i*_dimxy+j*_dimx+k]!=0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+j*_dimx+k, 1));
				continue;
			}


			//para of A(i+1,j,k) & A(i-1,j,k)
			//if z=0 boundary	
			if(i==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, (i+1)*_dimxy+j*_dimx+k, 1 + GetBoundaryPara(_boundaryType[4])));
			}
			//if z=dimz-1 boundary
			else if(i==_dimz-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, (i-1)*_dimxy+j*_dimx+k, 1 + GetBoundaryPara(_boundaryType[5])));
			}
			//no boundary in z-axis direction
			else{
				int iplus = 0;
				int iminus = 0;

				if(_holeMaskMap[(i+1)*_dimxy+j*_dimx+k]==0)//Normal
					iplus += 1;
				else if(_holeMaskMap[(i+1)*_dimxy+j*_dimx+k]==1)//Dirichlet
					;//iplus = 0;
				else if(_holeMaskMap[(i+1)*_dimxy+j*_dimx+k]==2)//Neumann
					iminus += 1;
				
				if(_holeMaskMap[(i-1)*_dimxy+j*_dimx+k]==0)//Normal
					iminus += 1;
				else if(_holeMaskMap[(i-1)*_dimxy+j*_dimx+k]==1)//Dirichlet
					;//iminus = 0;
				else if(_holeMaskMap[(i-1)*_dimxy+j*_dimx+k]==2)//Neumann
					iplus += 1;				

				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, (i+1)*_dimxy+j*_dimx+k, iplus));
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, (i-1)*_dimxy+j*_dimx+k, iminus));
			}
			
			//para of A(i,j+1,k) & A(i,j-1,k)
			//if y=0 boundary
			if(j==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+(j+1)*_dimx+k, 1 + GetBoundaryPara(_boundaryType[2]) ));
			}
			//if y=dimy-1 boundary
			else if(j==_dimy-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+(j-1)*_dimx+k, 1 + GetBoundaryPara(_boundaryType[3]) ));
			}
			//no boundary in y-axis direction
			else{
				int jplus = 0;
				int jminus = 0;

				if(_holeMaskMap[i*_dimxy+(j+1)*_dimx+k]==0)//Normal
					jplus += 1;
				else if(_holeMaskMap[i*_dimxy+(j+1)*_dimx+k]==1)//Dirichlet
					;//jplus = 0;
				else if(_holeMaskMap[i*_dimxy+(j+1)*_dimx+k]==2)//Neumann
					jminus += 1;
				
				if(_holeMaskMap[i*_dimxy+(j-1)*_dimx+k]==0)//Normal
					jminus += 1;
				else if(_holeMaskMap[i*_dimxy+(j-1)*_dimx+k]==1)//Dirichlet
					;//jminus = 0;
				else if(_holeMaskMap[i*_dimxy+(j-1)*_dimx+k]==2)//Neumann
					jplus += 1;				

				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+(j+1)*_dimx+k, jplus));
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+(j-1)*_dimx+k, jminus));
			}

			//para of A(i,j,k+1) & A(i,j,k-1)
			//if x=0 boundary
			if(k==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+j*_dimx+k+1, 1 + GetBoundaryPara(_boundaryType[0]) ));
			}
			//if x=dimx-1 boundary
			else if(k==_dimx-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+j*_dimx+k-1, 1 + GetBoundaryPara(_boundaryType[1]) ));
			}
			//no boundary in x-axis direction
			else{
				int kplus = 0;
				int kminus = 0;

				if(_holeMaskMap[i*_dimxy+j*_dimx+k+1]==0)//Normal
					kplus += 1;
				else if(_holeMaskMap[i*_dimxy+j*_dimx+k+1]==1)//Dirichlet
					;//kplus = 0;
				else if(_holeMaskMap[i*_dimxy+j*_dimx+k+1]==2)//Neumann
					kminus += 1;
				
				if(_holeMaskMap[i*_dimxy+j*_dimx+k-1]==0)//Normal
					kminus += 1;
				else if(_holeMaskMap[i*_dimxy+j*_dimx+k-1]==1)//Dirichlet
					;//kminus = 0;
				else if(_holeMaskMap[i*_dimxy+j*_dimx+k-1]==2)//Neumann
					kplus += 1;				

				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+j*_dimx+k+1, kplus));
				coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+j*_dimx+k-1, kminus));
			}
			
			//para of A(i,j,k)
			coefficients.push_back(Eigen::Triplet<double>( i*_dimxy+j*_dimx+k, i*_dimxy+j*_dimx+k, -6));//center));
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
	Eigen::VectorXd b(_dimx*_dimy*_dimz);	
	int _dimxy = _dimx * _dimy;

	//fill vector with force function
	for(int i=0;i<_dimz;i++)
		for(int j=0;j<_dimy;j++)
			for(int k=0;k<_dimx;k++)
			{	
				if(_holeMaskMap[i*_dimxy+j*_dimx+k] == 0)
					b(i*_dimxy+j*_dimx+k)=_forceFunc(i,j,k,_dimz,_dimy,_dimx)*_dx*_dx;
			}

	int i,j,k;
	//use different boundary condition

	//z=0 boundary
	i=0;
	for(j = 0; j < _dimy; j++)
	for(k = 0; k < _dimx; k++){
		if(_holeMaskMap[i*_dimxy+j*_dimx+k]==0)
		if(_boundaryType[4] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (-1,j,k)
			b(i*_dimxy+j*_dimx+k)-=_boundaryFunc[4](-1, j, k, _dimz, _dimy, _dimx);
		else if(_boundaryType[4] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (0,j,k)
			b(i*_dimxy+j*_dimx+k)+=_boundaryFunc[4](0, j, k, _dimz, _dimy, _dimx) * _dx * 2;
	}

	//z=dimz-1 boundary
	i=_dimz-1;
	for(j = 0; j < _dimy; j++)
	for(k = 0; k < _dimx; k++){
		if(_holeMaskMap[i*_dimxy+j*_dimx+k]==0)
		if(_boundaryType[5] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (dimz,j,k)
			b(i*_dimxy+j*_dimx+k)-=_boundaryFunc[5](_dimz, j, k, _dimz, _dimy, _dimx);
		else if(_boundaryType[5] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (dimz-1,j,k)
			b(i*_dimxy+j*_dimx+k)-=_boundaryFunc[5](_dimz-1, j, k, _dimz, _dimy, _dimx) * _dx * 2;
	}

	//y=0 boundary
	j=0;
	for(i = 0; i < _dimz; i++)
	for(k = 0; k < _dimx; k++){
		if(_holeMaskMap[i*_dimxy+j*_dimx+k]==0)
		if(_boundaryType[2] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (i,-1,k)
			b(i*_dimxy+j*_dimx+k)-=_boundaryFunc[2](i, -1, k, _dimz, _dimy, _dimx);
		else if(_boundaryType[2] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (i,0,k)
			b(i*_dimxy+j*_dimx+k)+=_boundaryFunc[2](i, 0, k, _dimz, _dimy, _dimx) * _dx * 2;
	}

	//y=dimy-1 boundary
	j=_dimy-1;
	for(i = 0; i < _dimz; i++)
	for(k = 0; k < _dimx; k++){
		if(_holeMaskMap[i*_dimxy+j*_dimx+k]==0)
		if(_boundaryType[3] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (i,dimy,k)
			b(i*_dimxy+j*_dimx+k)-=_boundaryFunc[3](i, _dimy, k, _dimz, _dimy, _dimx);
		else if(_boundaryType[3] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (i,dimy-1,k)
			b(i*_dimxy+j*_dimx+k)-=_boundaryFunc[3](i, _dimy-1, k, _dimz, _dimy, _dimx) * _dx * 2;
	}

	k=0;
	//x=0 boundary
	for(i = 0; i < _dimz; i++)
	for(j = 0; j < _dimy; j++){
		if(_holeMaskMap[i*_dimxy+j*_dimx+k]==0)
		if(_boundaryType[0] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (i,j,-1)
			b(i*_dimxy+j*_dimx+k)-=_boundaryFunc[0](i, j, -1, _dimz, _dimy, _dimx);
		else if(_boundaryType[0] & PB_BOUND_Neumann)
			//Neumann, boundary is the gradient on (i,j,0)
			b(i*_dimxy+j*_dimx+k)+=_boundaryFunc[0](i, j, 0, _dimz, _dimy, _dimx) * _dx * 2;
	}

	//x=dimx-1 boundary
	k=_dimx-1;
	for(i = 0; i < _dimz; i++)
	for(j = 0; j < _dimy; j++){
		if(_holeMaskMap[i*_dimxy+j*_dimx+k]==0)
		if(_boundaryType[1] & PB_BOUND_Dirichlet)
			//Dirichlet, boundary is the value on (i,j,dimx)
			b(i*_dimxy+j*_dimx+k)-=_boundaryFunc[1](i, j, _dimx, _dimz, _dimy, _dimx);
		else if(_boundaryType[1] & PB_BOUND_Neumann)
			//Neumann boundary is the gradient on (i,j,dimx-1)
			b(i*_dimxy+j*_dimx+k)+=_boundaryFunc[1](i, j, _dimx-1, _dimz, _dimy, _dimx) * _dx * 2;
	}


	for(i = 0; i < _dimz; i++)
	for(j = 0; j < _dimy; j++)
	for(k = 0; k < _dimx; k++)
			if(_holeMaskMap[i*_dimxy+j*_dimx+k] == 0){//not hole
				//z-1 neighbour is hole
				if(i-1>0 && _holeMaskMap[(i-1)*_dimxy+j*_dimx+k])

					if(_holeMaskMap[(i-1)*_dimxy+j*_dimx+k]==PB_BOUND_Dirichlet)
						b(i*_dimxy+j*_dimx+k)-=_holeBoundaryMap[(i-1)*_dimxy+j*_dimx+k];
					else if(_holeMaskMap[(i-1)*_dimxy+j*_dimx+k]==PB_BOUND_Neumann)
						b(i*_dimxy+j*_dimx+k)+= _dx * 2 * _gradientZFunc(i,j,k,_dimz,_dimy,_dimx);

				//z+1 neighbour is hole
				if(i+1<_dimz && _holeMaskMap[(i+1)*_dimxy+j*_dimx+k])

					if(_holeMaskMap[(i+1)*_dimxy+j*_dimx+k]==PB_BOUND_Dirichlet)
						b(i*_dimxy+j*_dimx+k)-=_holeBoundaryMap[(i+1)*_dimxy+j*_dimx+k];
					else if(_holeMaskMap[(i+1)*_dimxy+j*_dimx+k]==PB_BOUND_Neumann)
						b(i*_dimxy+j*_dimx+k)-= _dx * 2 * _gradientZFunc(i,j,k,_dimz,_dimy,_dimx);

				//y-1 neighbour is hole
				if(j-1>0 && _holeMaskMap[i*_dimxy+(j-1)*_dimx+k])

					if(_holeMaskMap[i*_dimxy+(j-1)*_dimx+k]==PB_BOUND_Dirichlet)
						b(i*_dimxy+j*_dimx+k)-=_holeBoundaryMap[i*_dimxy+(j-1)*_dimx+k];
					else if(_holeMaskMap[i*_dimxy+(j-1)*_dimx+k]==PB_BOUND_Neumann)
						b(i*_dimxy+j*_dimx+k)+= _dx * 2 * _gradientYFunc(i,j,k,_dimz,_dimy,_dimx);

				//y+1 neighbour is hole
				if(j+1<_dimy && _holeMaskMap[i*_dimxy+(j+1)*_dimx+k])

					if(_holeMaskMap[i*_dimxy+(j+1)*_dimx+k]==PB_BOUND_Dirichlet)
						b(i*_dimxy+j*_dimx+k)-=_holeBoundaryMap[i*_dimxy+(j+1)*_dimx+k];
					else if(_holeMaskMap[i*_dimxy+(j+1)*_dimx+k]==PB_BOUND_Neumann)
						b(i*_dimxy+j*_dimx+k)-= _dx * 2 * _gradientYFunc(i,j,k,_dimz,_dimy,_dimx);

				//x-1 neighbour is hole
				if(k-1>0 && _holeMaskMap[i*_dimxy+j*_dimx+k-1])

					if(_holeMaskMap[i*_dimxy+j*_dimx+k-1]==PB_BOUND_Dirichlet)
						b(i*_dimxy+j*_dimx+k)-=_holeBoundaryMap[i*_dimxy+j*_dimx+k-1];
					else if(_holeMaskMap[i*_dimxy+j*_dimx+k-1]==PB_BOUND_Neumann)
						b(i*_dimxy+j*_dimx+k)+= _dx * 2 * _gradientXFunc(i,j,k,_dimz,_dimy,_dimx);

				//x+1 neighbour is hole
				if(k+1<_dimx && _holeMaskMap[i*_dimxy+j*_dimx+k+1])

					if(_holeMaskMap[i*_dimxy+j*_dimx+k+1]==PB_BOUND_Dirichlet)
						b(i*_dimxy+j*_dimx+k)-=_holeBoundaryMap[i*_dimxy+j*_dimx+k+1];
					else if(_holeMaskMap[i*_dimxy+j*_dimx+k+1]==PB_BOUND_Neumann)
						b(i*_dimxy+j*_dimx+k)-= _dx * 2 * _gradientXFunc(i,j,k,_dimz,_dimy,_dimx);
			}else{//is a hole
				b(i*_dimxy+j*_dimx+k) = _holeBoundaryMap[i*_dimxy+j*_dimx+k];
			}

	return b;
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
double sinpowbc(int i, int total){
	Eigen::ArrayXd arr= Eigen::ArrayXd::LinSpaced(total, -SQRTPI,SQRTPI).pow(2);
	//printf("%f\t",arr(i));
	return sin(arr(i)+(SQRTPI*(total+1)/(total-1))*(SQRTPI*(total+1)/(total-1)));
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

//z = x+y+z
double plusof(int zi, int yi, int xi, int dimz, int dimy, int dimx){
	return xi+yi+zi;
}

double plusff(int zi, int yi, int xi, int dimz, int dimy, int dimx){
	return 0;
}

double plusbc(int zi, int yi, int xi, int dimz, int dimy, int dimx){
	return xi+yi+zi;
}

double plusgx(int zi, int yi, int xi, int dimz, int dimy, int dimx){
	return 1;
}
double plusgy(int zi, int yi, int xi, int dimz, int dimy, int dimx){
	return 1;
}
double plusgz(int zi, int yi, int xi, int dimz, int dimy, int dimx){
	return 1;
}