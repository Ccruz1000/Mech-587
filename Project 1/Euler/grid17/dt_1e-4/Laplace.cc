#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include "Base.h"

double f(double x, double y);
double f2(double x, double y);
void initializeUe(Vector &ue, const Grid &G);
void initializeU(Vector &u, const Grid &G);


/*================================================================================================
 * The current code can be compiled successfully without modification. You need to do this first. 
 * However, it gives you -nan in the calculation.
 * To make it right, you need to correct the following four functions. 
 * You can find them later in this file
 *================================================================================================
*/

void computeMatrix(Matrix &M, const Grid &G);
void computeTransientMatrix(Matrix &M, const Grid &G, const double &dt);
void computeDiffusion(Vector &R, const Matrix &M, const Vector &u);
void applyBC(Vector &R, Vector &du, const Grid &G);

/*===============================================================================================*/

// void storeVTKStructured(const Vector &b, const Grid &G, const char *fileName);

double f(double x, double y)
{
	return pow(x,4) + pow(y,4)- 6*pow(x,2)*pow(y,2);
}

double f2(double x, double y)
{
	return exp(-100*(pow(x,2) + pow(y,2)));
}

void initializeUe(Vector &ue, const Grid &G)
{
	unsigned long i,j;
	double x, y;

	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++){
			x = i*dx, y=j*dy;
			ue(i,j) = f(x,y);
		}
}

void initializeU(Vector &u, const Grid &G)
{
	unsigned long i,j;
	double x,y;

	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();

	for(i = 0; i < Nx; i++){
		x = i*dx;
		
		j = 0; y = j*dy;
		u(i,j) = f(x,y);

		j = Ny-1; y = j*dy;
		u(i,j) = f(x,y);
	}

	for(j = 0; j < Ny; j++){
		y = j*dy;

		i = 0; x = i*dx;
		u(i,j) = f(x,y);

		i = Nx-1; x = i*dx;
		u(i,j) = f(x,y);
	}

	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++){
			x = i*dx, y = j*dy;
			u(i,j) = f2(x,y);
		}
}

/*================================================================================
 * This is the function for constructing the LHS matrix for steady-state solution
 * The variables you need to correct are:
 * a, b, M(i,j,0), M(i,j,1), M(i,j,2), M(i,j,3), M(i,j,4)
 * You can use dx, dy and operators +,-,*,/ in the calculation
 *================================================================================
*/

void computeMatrix(Matrix &M, const Grid &G)
{
			/*
			 *=====================================================================
			 *Added boundary conditions based on whats given in steady
			 *=====================================================================
			*/

	unsigned long i,j;
	const double dx = G.dx();
	const double dy = G.dy();

	const double a = 1.0;
	const double b = 0.0;

	/* 
	 *a is 1.0 as M(i, j, 2) is the value that solves the point we are looking at. Along the
	 *boundary, this is 0. Seen in denotation for sparse matrix
	*/ 


	for(i = 1; i < G.Nx()-1; i++)
		for(j = 1; j < G.Ny()-1; j++){
			M(i,j,0) = 1.0 / (dx * dx);
			M(i,j,1) = 1.0 / ( dy * dy);
			M(i,j,2) = - 2.0 / (dx * dx) - 2.0 / (dy * dy);
			M(i,j,3) = 1.0 / (dy * dy);
			M(i,j,4) = 1.0 / (dx * dx);
		}

	for(i = 0; i < G.Nx(); i++)
		for(int t = 0; t < 5; t++){
			M(i,0,t) = (t == 2 ? a : b);  // Bottom boundary
			M(i,G.Ny()-1, t) = (t == 2 ? a : b);  // Top boundary
		}

	for(j = 0; j < G.Ny(); j++)
		for(int t = 0; t < 5; t++){
			M(0,j,t) = (t == 2 ? a : b);  // Left Boundary
			M(G.Nx()-1,j,t) = (t == 2 ? a : b);  // Right Boundary
		}
}

/*============================================================================
 * This is the function for constructing the LHS matrix for transient solution
 * You need to set up the values in two different ways: 
 * (1) Explicit Euler (2) Trapezoidal rule
 * The variables you need to correct are:
 * a, b, M(i,j,0), M(i,j,1), M(i,j,2), M(i,j,3), M(i,j,4)
 * You can use dx, dy, dt and operators +,-,*,/ in the calculation
 *============================================================================
*/

void computeTransientMatrix(Matrix &M, const Grid &G, const double &dt)
{			
			/*
			 *=====================================================================
			 *Added boundary conditions based on whats given in transient
			 *=====================================================================
			*/
	unsigned long i,j;
	const double dx = G.dx();
	const double dy = G.dy();
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	const double a = 1.0;
	const double b = 0.0;
	
	/* 
	 *a is 1.0 as M(i, j, 2) is the value that solves the point we are looking at. Along the
	 *boundary, this is 0. Seen in denotation for sparse matrix
	*/ 

	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++) {
			/*
			Boundary conditions the same here, but with the difference (I/dt)
			I is 1 along the diagnol, and 0 everywhere else. Therefore if the current position is not along the diagnol I/dt becomes 0, 
			and if it is, it becomes 1/dt. Only M(i, j, 2) is along the diagnol.
			*/
			M(i,j,0) = 0.0;
			M(i,j,1) = 0.0;
			M(i,j,2) = 1/dt;
			M(i,j,3) = 0.0;
			M(i,j,4) = 0.0;
		}

	for(i = 0; i < Nx; i++)
		for(int t = 0; t < 5; t++){
			M(i,0,t) = (t == 2 ? a : b);
			M(i,Ny-1, t) = (t == 2 ? a : b);
		}

	for(j = 0; j < Ny; j++)
		for(int t = 0; t < 5; t++){
			M(0,j,t) = (t == 2 ? a : b);
			M(Nx-1,j,t) = (t == 2 ? a : b);
		}
}

/*===================================================================================
 * This is the function for constructing the residual vector R=Au
 * The variables you need to correct is:
 * R(i,j)
 * You can use dx, dy, u(i-1,j), u(i,j-1), u(i,j) u(i,j+1) u(i+1,j)
 * and operators +,-,*,/ in the calculation
 *===================================================================================
*/

void computeDiffusion(Vector &R, const Vector &u, const Grid &G)
{
	// Added boundary condition calculated using Central Differencing method
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();

	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++)
			// Computed using central differencing method for both X and Y spatial derivative
			R(i,j) = (((u(i-1, j) - (2 * u(i, j)) + u(i+1, j)) / (dx * dx)) + ((u(i, j-1) - (2 * u(i, j)) + u(i, j+1)) / (dy * dy)));

}

/*===================================================================================
 * This is the function for imposing bounadry condition
 * The variables you need to correct is:
 * a
 *===================================================================================
*/
void applyBC(Vector &R, Vector &du, const Grid &G)
{
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();

	const double a = 0.0;
	// The above value is a as we have a Dirichlet boundary condition which remains as a constant value. Therefore we don't 
	// need to apply the boundary conditions again. 

	for(i = 0; i < Nx; i++){
		// Bottom Boundary
		j = 0;
		R(i,j) = du(i,j) = a;
		// Top boundary
		j = Ny-1;
		R(i,j) = du(i,j) = a;
	}

	for(j = 0; j < Ny; j++){
		// Left Boundary 
		i = 0;
		R(i,j) = du(i,j) = a;
		// Right boundary 
		j = Ny-1;
		R(i,j) = du(i,j) = a;
	}
}

const Vector solveSteadyLaplace(const Grid &G)
{
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	Matrix A(Nx,Ny);
	Vector R(Nx,Ny);
	Vector u(Nx,Ny);
	Vector du(Nx,Ny);
	
	//Initialize timing variables
	clock_t steady_start;
	double  steady_duration;

	//Initialize A,b,u,du,ue
	initializeU(u, G);
	computeMatrix(A, G);
	char matName[20] = "SteadyA.dat";
	A.storeA(matName);

	//Compute Residual and Residual norm
	double R0 = 0.0;
	computeDiffusion(R, u, G);
	applyBC(R,du,G);
	R = -R;
	R0 = R.L2Norm();
	char vecName[20] = "SteadyR0.dat";
	R.storeV(vecName);
	printf("Initial Residual Norm = %14.12e\n", R0);

	//Start Outer Loop
	int m = 0;

	steady_start = clock();

	while(m < 10){
		//Solve the linear system Adu = -R
		solveGS(du,A,R);	

		//Update the solution u
		u = u + du;

		//Compute the Residual
		computeDiffusion(R, u, G);
		applyBC(R, du, G);
		R = -R;
		double R1 = R.L2Norm();

		printf("Outer Iteration = %d\n",++m); 
		printf("Residual Norm = %14.12e,\n Residual Norm Ratio (R/R0) = %14.12e\n", R1, R1/R0);

		//Check convergence
		if(R1/R0 < 1e-5)
			break;
	}

	steady_duration = (clock() - steady_start) / (double) CLOCKS_PER_SEC;
	printf("Steady-state Solution time = %14.12e\n", steady_duration);
	
	char fname[20] = "steadyPhi.vtk";
	storeVTKStructured(u, G, fname);

	return u;
}

const Vector solveUnsteadyLaplace(const Grid &G, const double dt)
{
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	Matrix A(Nx,Ny);
	Vector R(Nx,Ny);
	Vector u(Nx,Ny);
	Vector du(Nx,Ny);
	
	//Initialize timing variables
	clock_t unsteady_start;
	double  unsteady_duration;

	//Initialize A,b,u,du,ue
	initializeU(u, G);
	computeTransientMatrix(A, G, dt);

	//Compute Residual and Residual norm
	double R0 = 0.0;
	computeDiffusion(R, u, G);
	applyBC(R,du,G);
	R0 = R.L2Norm();
	printf("Initial Residual Norm = %14.12e\n", R0);


	//Start Time Loop
	int itime = 0;
	double tf = 20;

	unsteady_start = clock();

	while(itime*dt <= tf){
		//Solve the linear system Adu = -R
		solveGS(du,A,R);	

		//Update the solution u
		u = u + du;

		//Compute the Residual
		computeDiffusion(R, u, G);
		applyBC(R, du, G);
		double R1 = R.L2Norm();

		printf("Time-Step = %d\n",++itime); 
		printf("Residual Norm = %14.12e,\n Residual Norm Ratio (R/R0) = %14.12e\n", R1, R1/R0);

		//Check convergence
		if(du.L2Norm() < 1e-8){
			printf("Steady state reached in %d time steps.\n Final time = %lf.\n",itime,itime*dt);
			break;
		}
	}

	unsteady_duration = (clock() - unsteady_start) / (double) CLOCKS_PER_SEC;
	printf("Transient Solution time = %14.12e\n", unsteady_duration);

	char fname[20] = "unsteadyPhi.vtk";
	storeVTKStructured(u, G, fname);

	return u;
}

void postProc(const Grid &G, const Vector &uT, const Vector &uS)
{
	unsigned long Nx, Ny, i, j;
	Nx = G.Nx();
	Ny = G.Ny();

	//Set up exact solution
	Vector ue(Nx,Ny);
	initializeUe(ue, G);
	FILE *postFile;
	postFile = fopen("TransSolError.dat","a");

	//Compute error for transient solution
	Vector e(Nx, Ny);
	e = ue - uT;
	printf("Transient Solution Error = %14.12e\n", e.L2Norm());
	fprintf(postFile,"%lu\t%14.12e\t%14.12e\n", G.size(), e.L2Norm(), e.LinfNorm(i,j));
	fclose(postFile);

	//Compute error for steady solution
	postFile = fopen("SteadySolError.dat","a");
	e = ue - uS;
	printf("Steady Solution Error = %14.12e\n", e.L2Norm());
	fprintf(postFile,"%lu\t%14.12e\t%14.12e\n", G.size(), e.L2Norm(), e.LinfNorm(i,j));
	fclose(postFile);

	//Compute difference error
	e = uT - uS;
	double Linf = e.LinfNorm(i,j); 
	printf("Difference Error (L2, Linf) = %14.12e, %14.12e at (%lu, %lu)\n", e.L2Norm(), Linf, i, j);

	char exactSolution[10] = "Ue.vtk";
	storeVTKStructured(ue, G, exactSolution);
}



int main()
{
	unsigned long Nx, Ny, i,j;
	double dt;
	
	FILE *inpFile;
	inpFile = fopen("Laplace.inp","r");

	fscanf(inpFile,"%lu %lu %lf", &Nx, &Ny, &dt);
	Grid G(Nx, Ny);
	
	const Vector uT = solveUnsteadyLaplace(G,dt);
	const Vector uS = solveSteadyLaplace(G);

	postProc(G, uT,uS);



}
