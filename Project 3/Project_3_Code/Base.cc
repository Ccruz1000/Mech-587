#include "Base.h"

void Operations::alloc1D(double **V, const unsigned long &N, bool isInit, double def) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}

	*V = new double[N];

	if(isInit)
		for(unsigned long i = 0; i < N; i++)
			(*V)[i] = 0.0;
}

void Operations::deAlloc1D(double **V, const unsigned long &N) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}
}

void Operations::copyArray(double *S, double *T, const unsigned long &N) const
{
	for(int i = 0; i < N; i++)
		T[i] = S[i];
}

void Operations::addVec(double *R, const double *V1, const double *V2, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = V1[i] + V2[i];
}

void Operations::subtractVec(double *R, const double *V1, const double *V2, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = V1[i] - V2[i];
}

void Operations::scaleVec(double *R, const double &a, const double *V1, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = a*V1[i];
}

void Operations::dotVec(double &R, const double *V1, const double *V2, const unsigned long &N) const
{
	R = 0.0;
	for(unsigned long i = 0; i < N; i++)
		R = R + (V1[i]*V2[i]);
}

void Operations::findAbsMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = -1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(fabs(V2[i]) > R){
			idx = i;
			R = fabs(V2[i]);
		}
}

void Operations::findMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = -1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(V2[i] > R){
			idx = i;
			R = V2[i];
		}
}

void Operations::findMin(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = 1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(V2[i] < R){
			idx = i;
			R = V2[i];
		}
}

Grid::Grid() : X(nullptr), Y(nullptr), nx(0), ny(0), N(0)
{}

Grid::Grid(unsigned long iNx, unsigned long iNy) : nx(iNx), ny(iNy), X(nullptr), Y(nullptr)
{
	xlim[0] = 0.0, xlim[1] = 1.0;
	ylim[0] = 0.0, ylim[1] = 1.0;
	N = nx*ny;
	dX = 1.0/(nx-1);
	dY = 1.0/(ny-1);

	O.alloc1D(&X, nx);
	O.alloc1D(&Y, ny);

	unsigned long i;
	for(i = 0; i < nx; i++)
		X[i] = i*dX;

	for(i = 0; i < ny; i++)
		Y[i] = i*dY;
}

Grid::Grid(unsigned long iNx, unsigned long iNy, double ixlim[2], double iylim[2]) : nx(iNx), ny(iNy), X(nullptr), Y(nullptr)
{
	xlim[0] = ixlim[0], xlim[1] = ixlim[1];
	ylim[0] = iylim[0], ylim[1] = iylim[1];
	N = nx*ny;
	dX = (xlim[1]-xlim[0])/(nx-1);
	dY = (ylim[1]-ylim[0])/(ny-1);

	O.alloc1D(&X, nx);
	O.alloc1D(&Y, ny);

	unsigned long i;
	for(i = 0; i < nx; i++)
		X[i] = xlim[0] + i*dX;

	for(i = 0; i < ny; i++)
		Y[i] = ylim[0] + i*dY;
}

Grid::~Grid()
{
	O.deAlloc1D(&X, nx);
	O.deAlloc1D(&Y, ny);
	nx = 0, ny = 0;
	N = 0;
	dX = 0.0, dY = 0.0;
}

void Grid::setGridSize(const unsigned long iNx, const unsigned long iNy)
{
	if(X != nullptr && Y != nullptr)
		return;

	nx = iNx, ny = iNy;
	dX = 1.0/(nx-1); dY = 1.0/(ny-1);
	N = nx*ny;
	O.alloc1D(&X,nx);
	O.alloc1D(&Y,ny);

	unsigned long i;
	for(i = 0; i < nx; i++)
		X[i] = i*dX;

	for(i = 0; i < ny; i++)
		Y[i] = i*dY;
}

Matrix::Matrix(): isInit(false)
{
	N = 0, Nx = 0, Ny = 0;
	for(int i = 0; i < 5; i++)
		A[i] = nullptr;
}

Matrix::Matrix(unsigned long iNx, unsigned long iNy) : Nx(iNx), Ny(iNy), N(iNx*iNy), isInit(false)
{
	for(int i = 0; i < 5; i++){
		A[i] = nullptr;
		O.alloc1D(&A[i], N);
	}
	isInit = true;

}

Matrix::~Matrix()
{
	for(int i = 0; i < 5; i++)
		if(A[i] != nullptr)
			O.deAlloc1D(&A[i], N);
	Nx = 0, Ny = 0, N = 0;
}

bool Matrix::isValid(unsigned short pos, unsigned long idxNode) const
{
	if(idxNode < N && (pos >= 0 && pos <= 4))
		return true;
	return false;
}

double& Matrix::operator()(unsigned long i, unsigned long j, unsigned short pos)
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	unsigned long idxNode = i*Ny+j;
	if(isValid(pos, idxNode))
		return A[pos][idxNode];
	else
	{
		printf("Attempt to access invalid memory. Please check. Exiting \n");
		exit(0);
		//Trigger exit code
	}
}

double Matrix::operator()(unsigned long i, unsigned long j, unsigned short pos) const
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check at (%lu,%lu). Exiting\n",i,j);
		exit(0);
	}

	unsigned long idxNode = i*Ny+j;
	if(isValid(pos, idxNode))
		return A[pos][idxNode];
	else
		return 0.0;
}

double& Matrix::operator()(unsigned long i, unsigned short pos)
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(pos, i))
		return A[pos][i];
	else
	{
		printf("Attempt to access invalid memory. Please check. Exiting \n");
		exit(0);
		//Trigger exit code
	}
}

double Matrix::operator()(unsigned long i, unsigned short pos) const
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(pos, i))
		return A[pos][i];
	else
		return 0.0;
}

void Matrix::setSize(unsigned long iNx, unsigned long iNy){
	Nx = iNx, Ny = iNy;
	N = Nx*Ny;
	for(int i = 0; i < 5; i++){
		A[i] = nullptr;
		O.alloc1D(&A[i], N);
	}
	isInit = true;
}

Vector::Vector():Nx(0), Ny(0), N(0), b(nullptr), isInit(false){}

Vector::Vector(unsigned long iNx, unsigned long iNy, bool initflag, double initVal): Nx(iNx), Ny(iNy), N(iNx*iNy), b(nullptr)
{
	O.alloc1D(&b,N,initflag,initVal);
	isInit = true;
}


Vector::Vector(const Vector &V1) : Nx(V1.Nx), Ny(V1.Ny), N(V1.N), isInit(false), b(nullptr)
{
	O.alloc1D(&b,N,false);
	isInit = true;
	O.copyArray(V1.b, b, N);
}

Vector::~Vector()
{
	O.deAlloc1D(&b,N);
	isInit = false;
	Nx = 0, Ny = 0; N = 0;
}

void Vector::setSize(unsigned long iNx, unsigned long iNy)
{
	Nx = iNx, Ny = iNy;
	N = Nx*Ny;
	b = nullptr;
	O.alloc1D(&b,N);
	isInit = true;
}

double& Vector::operator()(unsigned long i, unsigned long j)
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	unsigned long idxNode = i*Ny+j;
	if(isInit && isValid(idxNode))
		return b[idxNode];
	else{
		printf("attempted to access invalid memory location at (%lu,%lu). Please check. Exiting\n",i,j);
		exit(0);
	}
}

double Vector::operator()(unsigned long i, unsigned long j) const
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}
	unsigned long idxNode = i*Ny+j;
	if(isInit && isValid(idxNode))
		return b[idxNode];
	else
		return 0.0;
}


double& Vector::operator()(unsigned long i)
{
	if(!isInit){
		printf("Vector has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(i))
		return b[i];
	else{
		printf("attempted to access invalid memory location. Please check. Exiting\n");
		exit(0);
	}
}

double Vector::operator()(unsigned long i) const
{
	if(!isInit){
		printf("Vector has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(i))
		return b[i];
	else
		return 0.0;
}



Vector Vector::operator = (const Vector &V1)
{
	Nx = V1.Nx;
	Ny = V1.Ny;

	if(N == V1.N)
		O.copyArray(V1.b, b, V1.N);
	else{
		O.deAlloc1D(&b, N);
		N = V1.N;
		O.alloc1D(&b, N);
		O.copyArray(V1.b, b, N);
	}

	return *this;
}



Vector Vector::operator - ()
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b, -1.0, this->b, this->N);

	return R;
}

Vector Vector::operator + (Vector const &V)
{
	if(V.N != this->N){
		printf("Invalid vector addition operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny, false);
	O.addVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator - (const Vector &V)
{
	if(V.N != this->N){
		printf("Invalid vector subtraction operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny,false);
	O.subtractVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator * (const double &a)
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b,a,this->b,this->N);
	return R;
}

Vector Vector::operator + (Vector const &V) const
{
	if(V.N != this->N){
		printf("Invalid vector addition operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny, false);
	O.addVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator - (const Vector &V) const
{
	if(V.N != this->N){
		printf("Invalid vector subtraction operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny,false);
	O.subtractVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator * (const double &a) const
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b,a,this->b,this->N);
	return R;
}


Vector operator * (const double &a, const Vector &V1)
{
	Operations O;
	Vector R(V1.Nx, V1.Ny, false);
	O.scaleVec(R.b,a,V1.b,V1.N);
	return R;
}

const double operator % (const Vector &V1, const Vector &V2)
{
	double dotp = 0.0;
	Operations O;
	O.dotVec(dotp, V1.b, V2.b, V1.size());
	return dotp;
}

Vector operator * (const Matrix &A, const Vector &V)
{
	size_t Nx = V.Nx, Ny = V.Ny;
	Vector R(Nx, Ny, false);
	for(long int i = 0; i < A.size(); i++){
		R(i) = A(i,0) * V(i-Ny)
			 + A(i,1) * V(i-1)
			 + A(i,2) * V(i)
			 + A(i,3) * V(i+1)
			 + A(i,4) * V(i+Ny);
	}

	return R;
}

const double Vector::L2Norm()
{
	double L2;
	O.dotVec(L2,this->b, this->b, this->N);
	L2 = sqrt(L2/(this->N));
	return L2;
}

const double Vector::LinfNorm(unsigned long &ix, unsigned long &iy)
{
	double Linf = -1e20;
	unsigned long idx;
	O.findAbsMax(Linf, idx, this->b, this->N);
	ix = idx/Ny;
	iy = idx%Ny;
	return Linf;
}

const double Vector::GetMax()
{
	double max = -1e20;
	unsigned long idx;
	O.findMax(max,idx,this->b, this->N);
	return max;
}

const double Vector::GetMin()
{
	double min = 1e20;
	unsigned long idx;
	O.findMin(min, idx, this->b, this->N);
	return min;
}



void Vector::storeV(char filename[50])
{
	FILE *f;
	f = fopen(filename,"w");
	for(unsigned long i = 0; i < Nx; i++){
		for(unsigned long j = 0; j < Ny; j++)
			fprintf(f,"%lu %lu %10e\n",i, j, (*this)(i,j));
		fprintf(f,"\n");
	}
}

void Matrix::storeA(char filename[50])
{
	FILE *f;
	f = fopen(filename,"w");
	for(unsigned long i = 0; i < Nx; i++){
		for(unsigned long j = 0; j < Ny; j++)
			fprintf(f,"%lu %lu %10e %10e %10e %10e %10e \n",i, j, (*this)(i,j,0), (*this)(i,j,1), (*this)(i,j,2), (*this)(i,j,3), (*this)(i,j,4));
		fprintf(f,"\n");
	}
}

Solver::Solver():	solverCode(Config::Instance()->selectlinSolver()),
					solverPreConditioner(Config::Instance()->selectPreConditioner()),
					maxIterations(Config::Instance()->maxSolverIterations()),
					omegaSOR(Config::Instance()->orfSOR()),
					relResidual(Config::Instance()->solverRelResidual()),
					absResidual(Config::Instance()->solverAbsResidual()),
					alphaS(Config::Instance()->solverPreCondFactor()),
					Nx(Config::Instance()->Nx()),Ny(Config::Instance()->Ny()),
					R0(Nx,Ny), R(Nx,Ny), v(Nx,Ny), p(Nx,Ny), h(Nx,Ny), s(Nx,Ny),
					t(Nx,Ny), y(Nx,Ny), w(Nx,Ny), z(Nx,Ny),
					rho0(1.0), alpha(1.0), omega(1.0),
					LU(Nx,Ny)
{}	

Solver::~Solver() 
{}


void Solver::solveGS(Vector &u, const Matrix &A, const Vector &b)
{
	size_t Ny = b.sizeNy();
	unsigned long i;
	const Vector* const u0 = &u;

	int iter = 0;
	while(1){
		double L2Norm = 0.0;
		for(i = 0; i < A.size(); i++){
			double Res = 0.0;
			Res = Res + A(i,0) * (*u0)(i-Ny);
			Res = Res + A(i,1) * (*u0)(i-1);
			Res = Res + A(i,2) * (*u0)(i);
			Res = Res + A(i,3) * (*u0)(i+1);
			Res = Res + A(i,4) * (*u0)(i+Ny);

			Res = b(i) - Res;

			L2Norm += (Res*Res);
			u(i) = u(i) + omegaSOR*(1.0/A(i,2))*Res;
		}
		L2Norm = sqrt(L2Norm/b.size());
		++iter;
		
		if(L2Norm < absResidual || iter > maxIterations)
		{
			printf("Gauss Seidel Iterations converged\n");
			printf("iterations = %d, Error Residuals = %14.12e\n",iter,L2Norm);
			break;
		}
	}
}

void Solver::decomposeILU(const Matrix &A)
{
	size_t P;
	double alpha = solverPreConditioner == 2 ? alphaS : 1.0;
	double den1, den2;

	for(long int i= 0; i < Nx; i++)
		for(long int j = 0; j < Ny; j++){
			P = i*Ny+j;						

			den1 = alpha*(j-1 < 0 ? 0.0 : LU(i,j-1,3));
			den2 = alpha*(i-1 < 0 ? 0.0 : LU(i-1,j,4));

			LU(i,j,0) = A(i,j,0)/(1.0 + den1);				//Lw
			LU(i,j,1) = A(i,j,1)/(1.0 + den2);				//Ls
			LU(i,j,2) = A(i,j,2)
			          + LU(i,j,0)*den1
			          + LU(i,j,1)*den2
			          - LU(i,j,0)*(j-1 < 0 ? 0.0 : LU(i,j-1,4))
			          - LU(i,j,1)*(i-1 < 0 ? 0.0 : LU(i-1,j,3));	//Lp
			LU(i,j,3) = (A(i,j,3)-den1*LU(i,j,0))/LU(i,j,2);		//Un
			LU(i,j,4) = (A(i,j,4)-den2*LU(i,j,1))/LU(i,j,2);		//Ue
		}
}

void Solver::solveLU(Vector &u, const Matrix &LU, const Vector &b)
{
	//Forward Substitution - solve L(Ux) = b
	for(long int i = 0; i < Nx; i++)
		for(long int j = 0; j < Ny; j++)
			u(i,j) = (1.0/LU(i,j,2))
				   * (b(i,j) 
				   	- LU(i,j,0)*(i-1 < 0 ? 0.0 : u(i-1,j))
				   	- LU(i,j,1)*(j-1 < 0 ? 0.0 : u(i,j-1)));

	//Back Substitution - solve Ux = uf
	for(long int i = Nx-1; i >= 0; i--)
		for(long int j = Ny-1; j >= 0; j--)
			u(i,j) = u(i,j) 
				   - LU(i,j,3)*(j+1 > 0 ? 0.0 : u(i,j+1))
				   - LU(i,j,4)*(i+1 > 0 ? 0.0 : u(i+1,j));

}

void Solver::solveBiCG(Vector &u, const Matrix &A, const Vector &b)
{
	size_t Nx = Nx;
	size_t Ny = Ny;
	unsigned long iters = 0;

	//Initialize Residual
	R0 = b - A*u;

	R = R0;
	double R0Norm = R0.L2Norm(); double RNorm;
	printf("BiCGStab Solver Initial Norm = %14.12e\n",R0Norm);
	if(R0Norm < absResidual){
		printf("%lu %14.12e \n",iters, R0Norm);
		return;
	}

	//Setup parameters and vectors
	rho0 = alpha = omega = 1.0;

	if(solverPreConditioner == 0)
	{
		while(iters < maxIterations){
			RNorm = R.L2Norm();
			// Convergence condition for R
			if(RNorm < absResidual || RNorm/R0Norm < relResidual)
				break;

			rhoi = R0%R;
			beta = (rhoi/rho0)*(alpha/omega);
			rho0 = rhoi;

			p = R + beta*(p - omega*v);
			v = A*p;
			alpha = rhoi/(R0%v);
			h = u + alpha*p;
			//Convergence condition for h
			if(fabs(alpha)*p.L2Norm() < absResidual){
				u = h;
				break;
			}

			s = R - alpha*v;
			t = A*s;
			omega = (t%s)/(t%t);
			u = h + omega*s;
			// Convergence condition for u
			if(fabs(omega)*s.L2Norm() < absResidual)
				break;

			R = s - omega*t;
			iters++;
		}
		printf("BiCGStab Converged\n");
		printf("%lu %14.12e %14.12e\n",iters, RNorm, RNorm/R0Norm);
	}
	else if(solverPreConditioner == 1 || solverPreConditioner == 2)
	{
		//LU decomposition
		decomposeILU(A);
		printf("LU decomposition done\n");

		//Loop over the search process
		while(iters < maxIterations)
		{
			RNorm = R.L2Norm();
			// Convergence condition for R
			if(RNorm < absResidual || RNorm/R0Norm < relResidual)
				break;

			rhoi = R0%R;
			beta = (rhoi/rho0)*(alpha/omega);
			rho0 = rhoi;

			p = R + beta*(p - omega*v);
			solveLU(y,LU,p);
			v = A*y;
			alpha = rhoi/(R0%v);
			h = u + alpha*y;
			//Convergence condition for h
			if(fabs(alpha)*y.L2Norm() < absResidual){
				u = h;
				break;
			}

			s = R - alpha*v;
			t = A*s;
			solveLU(z,LU,s);
			t = A*z;
			solveLU(w,LU,t);
			omega = (w%z)/(w%w);
			u = h + omega*z;
			// Convergence condition for u
			if(fabs(omega)*z.L2Norm() < absResidual)
				break;

			R = s - omega*t;
			iters++;
		}
		printf("BiCGStab Converged\n");
		printf("%lu %14.12e %14.12e\n",iters, RNorm, RNorm/R0Norm);

	}
}

void Solver::solve(Vector &u, const Matrix &A, const Vector &b)
{
	if(solverCode == 1)
		solveGS(u,A,b);
	else if(solverCode == 2)
		solveBiCG(u,A,b);
}

void storeVTKStructured(const Vector &u, const Grid &G, const char *fileName) 
{
	FILE *vtkFile;
	vtkFile = fopen(fileName,"w");
	unsigned long i,j;
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	fprintf(vtkFile,"# vtk DataFile Version 2.0\n");
	fprintf(vtkFile,"TITLE = \"Quad data\"\n");
	fprintf(vtkFile,"ASCII\n");
	fprintf(vtkFile,"DATASET STRUCTURED_GRID\n");
	fprintf(vtkFile,"DIMENSIONS %lu %lu %d\n", Nx, Ny, 1);
	fprintf(vtkFile,"POINTS %lu FLOAT\n",Nx*Ny);
	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++){
			fprintf(vtkFile, "%lf %lf %lf\n",G.x(i), G.y(j), 0.0);
		}

	fprintf(vtkFile,"POINT_DATA %lu\n",Nx*Ny);
	fprintf(vtkFile,"SCALARS phi FLOAT 1\n");
	fprintf(vtkFile,"LOOKUP_TABLE default\n");
	for(i= 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
			fprintf(vtkFile,"%14.12lf\n",u(i,j));

	fclose(vtkFile);
}

void storeVTKSolution(const Vector &u, const Vector &v, const Vector &p, const Grid &G, const char *fileName) 
{;
	FILE *vtkFile;
	vtkFile = fopen(fileName,"w");
	unsigned long i,j;
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	fprintf(vtkFile,"# vtk DataFile Version 2.0\n");
	fprintf(vtkFile,"TITLE = \"Quad data\"\n");
	fprintf(vtkFile,"ASCII\n");
	fprintf(vtkFile,"DATASET STRUCTURED_GRID\n");
	fprintf(vtkFile,"DIMENSIONS %lu %lu %d\n", Nx, Ny, 1);
	fprintf(vtkFile,"POINTS %lu FLOAT\n",Nx*Ny);
	for(j = 0; j < Ny; j++)
		for(i = 0; i < Nx; i++){
			fprintf(vtkFile, "%f %f %f\n",G.x(i), G.y(j), 0.0);
		}

	fprintf(vtkFile,"POINT_DATA %lu\n",Nx*Ny);
	fprintf(vtkFile,"SCALARS pressure FLOAT 1\n");
	fprintf(vtkFile,"LOOKUP_TABLE default\n");
	for(j= 0; j < Ny; j++)
		for(i = 0; i < Nx; i++)
			fprintf(vtkFile,"%lf\n",p(i,j));

	fprintf(vtkFile,"VECTORS velocity FLOAT\n");
	for(j= 0; j < Ny; j++)
		for(i = 0; i < Nx; i++)
			fprintf(vtkFile,"%lf %lf %lf\n",u(i,j),v(i,j),0.0);

	fclose(vtkFile);
}

Config* Config::instance = nullptr;

Config::Config()
{
	char ch[50];

	//Open input file to read parameters
	inputfile = fopen("NS.inp","r");

	//Read Grid Size
	fscanf(inputfile,"%s%lu",ch,&nx);
	fscanf(inputfile,"%s%lu",ch,&ny);

	//Read N-S parameters
	fscanf(inputfile,"%s%lf",ch,&Dt);
	fscanf(inputfile,"%s%lf",ch,&re);
	fscanf(inputfile,"%s%lu",ch,&maxTimeSteps);
	fscanf(inputfile,"%s%lf",ch,&alphaP);
	fscanf(inputfile,"%s%lf",ch,&alphaU);
	fscanf(inputfile,"%s%lf",ch,&convergence);
	fscanf(inputfile,"%s%lu",ch,&maxOuterIterations);

	//Read linear Solver parameters
	fscanf(inputfile,"%s%hu",ch,&linSolver);
	fscanf(inputfile,"%s%lf",ch,&omega);
	fscanf(inputfile,"%s%lf",ch,&relResidual); relResidual = pow(10,-relResidual);
	fscanf(inputfile,"%s%lf",ch,&absResidual); absResidual = pow(10,-absResidual);
	fscanf(inputfile,"%s%lu",ch,&maxIterations);

	fscanf(inputfile,"%s%hu",ch,&preConditioner);
	fscanf(inputfile,"%s%lf",ch,&alphaS);

	fclose(inputfile);
	inputfile = nullptr;

	printf("Input file read successfully\n");
}
