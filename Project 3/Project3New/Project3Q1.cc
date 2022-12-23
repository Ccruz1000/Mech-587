#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include "Base.h"

#include <iostream>
// Input Variables
unsigned long Nx = 17;
unsigned long Ny = 17;
double conv_u = 1.0;  //Convective u velocity cx
double conv_v = 0.05;  //Convective v velocit cv
double Re = 10;
double dt = 0.01;
unsigned long nTimeSteps = 10000;

// Function for the initial value from project 2
double f1(double x, double y)
{
	return exp(-1500*(pow((x-0.5),2) + pow((y-0.5),2)));
}

// Function for the exact solution 
double fexact(double x, double y, double u, double v, double alpha)
{
	return (5.0*((1-exp(u*x/alpha))/(1-exp(u/alpha)))) + (0.1*((1-exp(y*v/alpha))/(1-exp(v/alpha))));
}



// Function to calculate the exact solution from project 2 
// Same as InitializePhiExact
void initVelexact(Vector &vel_exact, const Grid &G, const Vector &uconv, const Vector &vconv)
{
	unsigned long i, j;
	double x, y;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();
	double alpha = 1 / Re;
	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			x = G.x(i);
            y = G.y(j);
			vel_exact(i,j) = fexact(x, y, uconv(i, j), vconv(i, j), alpha);
			
		}
	char fname[50] = "vel_exact.vtk";
	storeVTKStructured(vel_exact, G, fname);
}

// Initialize the convective velocity Cx and Cy for part 1
// Same as InitializeVelCD
void InitVelConv(Vector &u, Vector &v, const Vector &p, const Grid &G)
{
	unsigned long i, j;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();


	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			 u(i, j) = conv_u;
			 v(i, j) = conv_v;
		}
	 char xname[50] = "Initx_Conv.vtk";
	 storeVTKStructured(u, G, xname);
	 char vname[50] = "Initv_Conv.vtk";
	 storeVTKStructured(v, G, vname);
}

// Here we are going to use the velocity vector as phi from project 2. This is because
// in the momentum equation we are trying to solve for the velocity
// Same as InitializePhiCD
 void InitVel(Vector &vel, const Grid &G, const Vector &uconv, const Vector &vconv)
 {
 	unsigned long i, j;
 	double x, y;
 	const unsigned long Nx = G.Nx();
 	const unsigned long Ny = G.Ny();
 	const double dx = G.dx();
 	const double dy = G.dy();
 	double alpha = 1 / Re;
	
 	for(i = 0; i < Nx; i++)
 		for(j = 0; j < Ny; j++)
 		{
 			x = G.x(i);
             y = G.y(j); 
 			if(i==0)
 			{
 				vel(i,j) = 0.1*((1-exp((vconv(i, j)*y)/alpha))/((1-exp(vconv(i, j)/alpha))));
 			}
 			else if(i==Nx-1)
 			{
 				vel(i,j) =5+0.1*((1-exp((vconv(i, j)*y)/alpha))/((1-exp(vconv(i, j)/alpha))));
 			}
             else if(j==0)
             {
                 x = G.x(i);
                 y = G.y(j+1); 
                 vel(i,j) = f1(x,y);
             }
             else if(j==Ny-1)
             {
                 x = G.x(i);
                 y = G.y(j-1); 
                 vel(i,j) = f1(x,y);
             }
 			else
 			{
 				vel(i,j) = f1(x,y);
 			}
 		}

 }

// Calculate first derivative using central difference method 
//DS2 Project 2
 void CenDiff_1st(Vector &Flux_Conv_Curr, const Vector &vel, const Vector &uconv, const Vector &vconv, const Grid &G)
 {
 	unsigned long i, j;
 	double x, y;

 	const unsigned long Nx = G.Nx();
 	const unsigned long Ny = G.Ny();
 	const double dx = G.dx();
 	const double dy = G.dy();

 	for(i = 1; i < Nx-1; i++)
 		for(j = 1; j < Ny-1; j++)

         {
 			Flux_Conv_Curr(i,j)=(uconv(i,j)/(2*dx))*(vel(i+1,j)-vel(i-1,j)) + (vconv(i,j)/(2*dy))*(vel(i,j+1)-vel(i,j-1));
 		}
 }

// Compute time derivative using Adams-Bashforth
 void abs2Exp(Vector &vel, const Vector &Flux_Conv_Curr, const Vector &Flux_Conv_Prev)
 {
 	vel = vel - ((1.5 * Flux_Conv_Curr * dt) - (0.5 * Flux_Conv_Prev * dt));
 }

// Calculate diffusive derivatives using central difference method
// computeDiffusion Project 2
 void CenDiffuse(Vector &diffuse, const Vector &vel, const Grid &G)
 {
 	unsigned long i,j;
 	unsigned long Nx = G.Nx();
 	unsigned long Ny = G.Ny();
 	double dx = G.dx();
 	double dy = G.dy();
     double alpha = 1 / Re;

 	for(i = 1; i < Nx-1; i++)
 		for(j = 1; j < Ny-1; j++)
 			diffuse(i,j)=alpha*(((vel(i-1,j)-2*vel(i,j)+vel(i+1,j))/(dx*dx))+((vel(i,j-1)-2*vel(i,j)+vel(i,j+1))/(dy*dy))); 
 }

 // Compute the transient matrix based on project 2
 void computeTransientMatrix(Matrix &M, const Grid &G)
 {
 	unsigned long i,j;
 	const double dx = G.dx();
 	const double dy = G.dy();
 	unsigned long Nx, Ny;
 	Nx = G.Nx(), Ny = G.Ny();

 	const double a = 1;
 	const double b = 0;

 	for(i = 1; i < Nx-1; i++)
 		for(j = 1; j < Ny-1; j++) {
 			M(i,j,0) = -1/(2*dx*dx);
 			M(i,j,1) = -1/(2*dy*dy);
 			M(i,j,2) = (1/dt)+(1/(dx*dx))+(1/(dy*dy));
 			M(i,j,3) = -1/(2*dy*dy);
 			M(i,j,4) = -1/(2*dx*dx);
 		}

 	for(i = 0; i < Nx; i++)  //NEUMANN BC
     {
             M(i,0,0) = 0;
 			M(i,0,1) = 0;
 			M(i,0,2) = -(1/(dy));
 			M(i,0,3) = (1/(dy));
 			M(i,0,4) = 0;

             M(i,Ny-1,0) = 0;
 			M(i,Ny-1,1) = -(1/(dy));
 			M(i,Ny-1,2) = (1/(dy));
 			M(i,Ny-1,3) = 0;
 			M(i,Ny-1,4) = 0;
 		}
    
 	for(j = 0; j < Ny; j++) //DIRICH BC
 		for(int t = 0; t < 5; t++){
 			M(0,j,t) = (t == 2 ? a : b);
 			M(Nx-1,j,t) = (t == 2 ? a : b);
 		}
 }


 void computeResidual(Vector &R, const Vector &Flux_Conv_Curr, const Vector &Flux_Conv_Prev, const Vector &diffuse)
 {

     R = - (1.5 * Flux_Conv_Curr) + (0.5 * Flux_Conv_Prev) + diffuse;

 }

// Apply the boundary conditions. Currently using the same boundary conditions as project 2
 void applyBC(Vector &R, Vector &dvel, const Grid &G, const Vector &vel, const Vector &uconv, const Vector &vconv)
 {
 	unsigned long i,j;
 	unsigned long Nx = G.Nx();
 	unsigned long Ny = G.Ny();

     double u = 1;
     double v = 0.05;
 	double alpha = 1 / Re; 

 	const double a = 0.0;

     double dx = G.dx();
 	double dy = G.dy();

     double h1; 
     double h2;

 	for(i = 0; i < Nx; i++){
 		j = 0;
 		h1 = 0.1*((-vconv(i, j))/(alpha*(1-exp(vconv(i, j)/alpha)))); 
 		h2 = 0.1*((-vconv(i, j)*exp(vconv(i, j)/alpha))/(alpha*(1-exp(vconv(i, j)/alpha))));
         R(i,j) = h1 - (vel(i,j+1) - vel(i,j))/dy;

 		j = Ny-1;
         R(i,j) = h2 - (vel(i,j) - vel(i,j-1))/dy;
 	}

 	for(j = 0; j < Ny; j++){
 		i = 0;
 		R(i,j) = dvel(i,j) = a;

 		i = Nx-1;
 		R(i,j) = dvel(i,j) = a;
 	}
 }

void SolveConvectionDiffusion(const Grid &G)
{
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();

	// Initialize vectors to store solutions
	Vector uconv(Nx, Ny); // Convective u vector
	Vector vconv(Nx, Ny); // Convective v vector
	Vector vel_exact(Nx, Ny); // Exact solution
	Vector u_sol(Nx, Ny); // u we are solving for 
	Vector v_sol(Nx, Ny); // v we are solving for 
	Vector p_sol(Nx, Ny); // Pressure we are solving for 
	InitVelConv(uconv, vconv, p_sol, G); // Initialize convective vector
	initVelexact(vel_exact, G, uconv, vconv); // Initialize exact solution vector

	// Initialize matrices
	Matrix Ax(Nx, Ny); // A matrix for X momentum
	//Matrix Ay(Nx, Ny); // A Matrix for Y momentum

	// Initialize residuals and derivative vectors for solver
	// X-momentum
	Vector Rx(Nx, Ny); // X-momentum residual
	Vector Flux_Conv_Currx(Nx, Ny); // Current X-momentum Flux
	Vector Flux_Conv_Prevx(Nx, Ny); // Previous X-momentum flux
	Vector diffusex(Nx, Ny); // Diffusion flux x-momentum
	Vector dvelx(Nx, Ny); // Increment dvel in X
	//Vector Ry(Nx, Ny); // Y-momentum residual
	//Vector Flux_Conv_Curry(Nx, Ny); // Current Y-momentum Flux
	//Vector Flux_Conv_Prevy(Nx, Ny); // Previous Y momentum flux
	//Vector diffusey(Nx, Ny); // Diffusion flux y-momentum

	// Initialize transient matrices
	computeTransientMatrix(Ax, G);
	//computeTransientMatrix(Ay, G);

	// Compute residuals and residual norms
	double R0x = 0.0; // initial x residual
	//double R0y = 0.0; // initial y residual
	// Solve x momentum
	CenDiffuse(diffusex, u_sol, G);
	CenDiff_1st(Flux_Conv_Currx, u_sol, uconv, vconv, G);
	Flux_Conv_Prevx = Flux_Conv_Currx;
	computeResidual(Rx, Flux_Conv_Currx, Flux_Conv_Prevx, diffusex);
	applyBC(Rx, dvelx, G, u_sol, uconv, vconv);
	R0x = Rx.L2Norm();
	printf("Initial Residual Norm = %14.12e\n", R0x);
	//CenDiffuse(diffusey, v_sol, G);
	unsigned long itime = 0; double time = 0.0; bool last = false;
	// Store initial vtk file
	char fileName[50] = "solutionx_0.vtk";
	storeVTKStructured(u_sol, G, fileName);
	while(itime < nTimeSteps)
	{
		Flux_Conv_Prevx = Flux_Conv_Currx;
		solveGS(dvelx, Ax, Rx);
		u_sol = u_sol + dvelx;

		// Calculate diffusion vector 
		CenDiffuse(diffusex, u_sol, G);
		// Calculate current convection vector;
		CenDiff_1st(Flux_Conv_Currx, u_sol, uconv, vconv, G);
		// Calculate Residual
		computeResidual(Rx, Flux_Conv_Currx, Flux_Conv_Prevx, diffusex);
		applyBC(Rx, dvelx, G, u_sol, uconv, vconv);

		double R1x = Rx.L2Norm();

		time = (last == true ? nTimeSteps:(itime + 1)*dt);
		if(itime % 1 == 0 || last)
		{
			printf("Time = %lf, maxvel = %14.12e, minvel = %14.12e\n",time,u_sol.GetMax(),u_sol.GetMin());
            printf("Residual Norm = %14.12e,\nResidual Norm Ratio (R/R0) = %14.12e\n", R1x, R1x/R0x);
			printf("Convergence criteria = %14.12e\n", dvelx.L2Norm());
            std::cout<<std::endl;
		}

		if(last)
		{
			break;
		}
		if((itime + 1) % 10 == 0)
		{
			sprintf(fileName, "solutionx_%lu.vtk", itime + 1);
			storeVTKStructured(u_sol, G, fileName);
		}
		//Check convergence
		if(dvelx.L2Norm() < 1e-8){
			printf("Steady state reached in %lu time steps.\n Final time = %lf.\n",itime,itime*dt);
			break;
		}

		++itime;
		if(itime*dt > nTimeSteps){
			dt = nTimeSteps - (itime-1)*dt;
			last = true;
		}
	}



}

int main()
{
	double xlim[2] = {0, 1};
	double ylim[2] = {0, 1};
	Grid G(Nx, Ny, xlim, ylim);
	SolveConvectionDiffusion(G);

	return 0;
}