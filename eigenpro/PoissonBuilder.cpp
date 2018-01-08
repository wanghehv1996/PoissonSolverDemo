#include "PoissonBuilder.h"


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

/******************************************************
 * Function:SetBoundary(int boundarytype, double (*boundaryFunc)(int,int))
 * 
 * Class: PoissonBuilder
 * Set boundary type & function.
 * 
 ******************************************************/
int PoissonBuilder::SetBoundary(int boundarytype, double (*boundaryFunc)(int,int)){
	//left boundary
	//store in _boundary[0]
	if(boundarytype && PB_BOUND_LEFT){
		_boundaryType[0] = boundarytype & 0x3;
		_boundaryFunc[0] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where x = 0 \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//right boundary
	//store in _boundary[1]
	if(boundarytype && PB_BOUND_RIGHT){
		_boundaryType[1] = boundarytype & 0x3;
		_boundaryFunc[1] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where x = colnum \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//up boundary
	//store in _boundary[2]
	if(boundarytype && PB_BOUND_UP){
		_boundaryType[2] = boundarytype & 0x3;
		_boundaryFunc[2] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where y = 0 \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
	//down boundary
	//store in _boundary[3]
	if(boundarytype && PB_BOUND_DOWN){
		_boundaryType[3] = boundarytype & 0x3;
		_boundaryFunc[3] = boundaryFunc;
		if(_debugInfo && GetDebugBoundaryTypeStr(boundarytype)){
			printf("PoissonBuilder: Add %s boundary where y = rownum \n", 
				GetDebugBoundaryTypeStr(boundarytype));
		}
	}
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
			//para of A(i+1,j) & A(i-1,j)
			//if up boundary	
			if(i==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i+1)*_dimx+j, 1 + GetBoundaryPara(_boundaryType[2])));
			}
			//if down boundary
			else if(i==_dimy-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i-1)*_dimx+j, 1 + GetBoundaryPara(_boundaryType[3])));
			}
			//no boundary in y-axis direction
			else{
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i+1)*_dimx+j, 1));
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, (i-1)*_dimx+j, 1));
			}
			
			//para of A(i,j+1) & A(i,j-1)
			//if left boundary
			if(j==0){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j+1, 1 + GetBoundaryPara(_boundaryType[0]) ));
			}
			//if right boundary
			else if(j==_dimx-1){
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j-1, 1 + GetBoundaryPara(_boundaryType[1]) ));
			}
			//no boundary in x-axis direction
			else{

				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j+1, 1));
				coefficients.push_back(Eigen::Triplet<double>( i*_dimx+j, i*_dimx+j-1, 1));
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
			b(i*_dimx+j)=_forceFunc(i,j,_dimy,_dimx);
		}
	}

	//use different boundary condition

	int i,j;

	j=0;
	//left boundary
	for(i = 0; i < _dimy; i++){
		if(_boundaryType[0] & PB_BOUND_Dirichlet)
			//Dirichlet
			b(i*_dimx+j)-=_boundaryFunc[0](i, _dimy);
		else if(_boundaryType[0] & PB_BOUND_Neumann)
			//Neumann
			b(i*_dimx+j)+=_boundaryFunc[0](i, _dimy) * _dx * 2;
	}

	//right boundary
	j=_dimx-1;
	for(i = 0; i < _dimy; i++){
		if(_boundaryType[1] & PB_BOUND_Dirichlet)
			//Dirichlet
			b(i*_dimx+j)-=_boundaryFunc[1](i, _dimy);
		else if(_boundaryType[1] & PB_BOUND_Neumann)
			//Neumann
			b(i*_dimx+j)+=_boundaryFunc[1](i, _dimy) * _dx * 2;
	}

	//up boundary
	i=0;
	for(j = 0; j < _dimx; j++){
		if(_boundaryType[2] & PB_BOUND_Dirichlet)
			//Dirichlet
			b(i*_dimx+j)-=_boundaryFunc[2](j, _dimx);
		else if(_boundaryType[2] & PB_BOUND_Neumann)
			//Neumann
			b(i*_dimx+j)+=_boundaryFunc[2](j, _dimx) * _dx * 2;
	}

	//down boundary
	i=_dimy-1;
	for(j = 0; j < _dimx; j++){
		if(_boundaryType[3] & PB_BOUND_Dirichlet)
			//Dirichlet
			b(i*_dimx+j)-=_boundaryFunc[3](j, _dimx);
		else if(_boundaryType[3] & PB_BOUND_Neumann)
			//Neumann
			b(i*_dimx+j)+=_boundaryFunc[3](j, _dimx) * _dx * 2;
	}

	return b;
}



//zero force funciton
double zeroff(int, int, int, int){
	return 0;
}

//sqrt(Pi)
#define SQRTPI 1.77245385091

//force function of 
//z = sin(x*x+y*y)
//laplacian(z) = 4*cos(x*x+y*y) - 4*sin(x*x+y*y)*(x*x+y*y)
double sinpowff(int xi, int yi, int dimx, int dimy){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx+2, -SQRTPI,SQRTPI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy+2, -SQRTPI,SQRTPI);
	double x = arrx(xi+1);
	double y = arry(yi+1);
	return 2*2*cos(x*x + y*y) - sin(x*x + y*y)*4*(x*x+y*y);

}

//boundary condition of
//z = sin(x*x+y*y)
double sinpowbc(int i, int total){
	Eigen::ArrayXd arr= Eigen::ArrayXd::LinSpaced(total, -SQRTPI,SQRTPI).pow(2);
	//printf("%f\t",arr(i));
	return sin(arr(i)+SQRTPI*SQRTPI);
}


//force function of 
//z = sin(x+y)
//laplacian(z) = -sin(x+y)
double sinff(int xi, int yi, int dimx, int dimy){
	Eigen::ArrayXd arrx= Eigen::ArrayXd::LinSpaced(dimx, 0,M_PI);
	Eigen::ArrayXd arry= Eigen::ArrayXd::LinSpaced(dimy, 0,M_PI);
	double x = arrx(xi);
	double y = arry(yi);
	return -sin(x+y);	
}

//boundary condition of
//z = sin(x+y)
double sinbc(int i, int total){
	Eigen::ArrayXd arr= Eigen::ArrayXd::LinSpaced(total, 0,M_PI);
	return sin(arr(i));
}