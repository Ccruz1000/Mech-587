#include "Base.h"
#include <cstring>
#include <string>
//TODO - Update initializephicd to include the boundary conditions properly 
//(The equation given does not include 0 and 1 so use boundary condition)
//TODO Neuman Boundary conditions
//


/*================================================================================================
 * The current code can be compiled successfully without modification with Base.h and Base.cc.
 * However, it is not solving anything
 * To make it right, you need to complete all the functions as requried below. 
 * You can find them later in this file
 * You are free to reuse any of the code in Project #1
 *================================================================================================
*/

/*================================================================================================
 * Correct the initialization function for the 2D convection equation f0
 *================================================================================================
 */
double f0(double x, double y);
/*================================================================================================*/

void InitializePhi(Vector &phi, const Grid &G);
void InitializeVel(Vector &u, Vector &v, const Grid &G);
/*================================================================================================
 * Correct the initialization function for the 2D convection diffusion equation f1 
 * pay attention to boundary values (use the value from the boundary condition)
 *================================================================================================
 */
double f1(double x, double y);
/*================================================================================================*/

void InitializePhiCD(Vector &phi, const Vector &u, const Vector &v, const Grid &G, double alpha);
void InitializeVelCD(Vector &u, Vector &v, const Grid &G);
void InitializePhiCD_Exact(const Grid &G, double alpha);

/*================================================================================================
 * Complete the function for applying convective operator
 * (1) first order upwind (FOU)
 * (2) central difference (CDS2)
 * (3) second order upwind (SOU)
 *================================================================================================
 */
void FOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G);
void CDS2(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G);
void SOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G);
/*================================================================================================*/


/*================================================================================================
 * Complete the function for time-integration
 * (1) forward Euler (eulerExp)
 * (2) Adams-Bashforth second order (ab2Exp)
 *================================================================================================
 */
void eulerExp(Vector &phi, const Vector &fc_Curr, double &dt);
void abs2Exp(Vector &phi, const Vector &fc_Curr, const Vector &fc_Prev, double &dt);
/*================================================================================================*/

// Declare functions from project 1 code 
void computeTransientMatrix(Matrix &M, const Grid &G, const double &dt);
void computeDiffusion(Vector &R, const Vector &u, const Grid &G);
void applyBC(Vector &R, Vector &dphi, const Vector &phi, const Vector &v, const Grid &G, double alpha);
void initializePhie(Vector &phie, const Grid &G);


void SolveConvection(const Grid &G, const double tf, double dt, const unsigned short scheme);

/*================================================================================================
 * Complete the function for solving convection diffusion equation
 * Variables which might be useful for you has been define inside the function
 *================================================================================================
 */
void SolveConvectionDiffusion(const Grid &G, double tf, double dt, const unsigned short scheme);
/*================================================================================================*/



double f0(double x, double y){
	double x_pow = (x - 0.25) * (x - 0.25);
	double y_pow = y * y;
	return 5.0 * exp(-1500.0 * (x_pow + y_pow));
}

void InitializePhi(Vector &phi, const Grid &G)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i,j;
	double x,y;

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
			phi(i,j) = f0(G.x(i), G.y(j));
}

void InitializeVel(Vector &u, Vector &v, const Grid &G)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	unsigned long i, j;
	double x,y;

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++){
			u(i,j) = G.y(j);
			v(i,j) = -G.x(i);
		}
}

double f1(double x, double y)
{	
	double xpow = x - 0.5;
	double ypow = y - 0.5;
	return 1.0 * exp(-1500 * ((xpow * xpow) + (ypow * ypow)));
}

void InitializePhiCD(Vector &phi, const Vector &u, const Vector &v, const Grid &G, double alpha)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i,j;

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			// Apply given boundary condition to left and right boundary 
			if (i == 0)
			{
				phi(i, j) = 0.1 * ((1 - exp(G.y(j) * v(i, j) / alpha)) / (1 - exp(v(i, j) / alpha)));
		
			}
			else if (i == (Nx - 1))
			{
				phi(i, j) = 5.0 + 0.1 * ((1 - exp(G.y(j) * v(i, j) / alpha)) / (1 - exp(v(i, j) / alpha)));
			}
			// Apply exact solution to top and bottom boundary 
			// else if(j == 0 || j == (Ny - 1))
			// {
			// 	phi(i, j) = 5.0 * ((1 - exp(G.x(i) * u(i, j) / alpha)) / (1 - exp(u(i, j) / alpha))) + 
			// 	            0.1 * ((1 - exp(G.y(j) * v(i, j) / alpha)) / (1 - exp(v(i, j) / alpha)));
		
			// }
			// Apply initial condition to internal points
		    else
		    {
				phi(i,j) = f1(G.x(i) , G.y(j));
			}
			
		}
}

void InitializeVelCD(Vector &u, Vector &v, const Grid &G)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	unsigned long i, j;
	double x,y;

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			u(i,j) = 1.0;
			v(i,j) = 0.05;
		}
}

void InitializePhiCD_Exact(const Grid &G, double alpha)
{
	// Initialize values for for loop
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i, j;
	
	// Initialize velocity and solution vectors
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);
	Vector phi(Nx, Ny);

	InitializePhiCD(phi, u, v, G, alpha);
	InitializeVelCD(u, v, G);

	// Solve exact equation
	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			phi(i, j) = 5.0 * (1 - exp(u(i, j) * G.x(i) / alpha)) / (1 - exp(u(i, j) / alpha)) + 
			0.1 * (1 - exp(v(i, j) * G.y(j) / alpha)) / (1 - exp(v(i, j) / alpha));
			
		}
	// Save VTK file
	char fname[20] = "PhiCD_Exact.vtk";
	storeVTKStructured(phi, G, fname);
}

void FOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
	unsigned long i, j;
	const double dx = G.dx();
	const double dy = G.dy();

    // Describe interior nodes
	for(i = 1; i < G.Nx() - 1; i++)
		for(j = 1; j < G.Ny() - 1; j++)
		{
			// This one doesnt work because abs returns an int
			/*fc_Curr(i, j) = (u(i, j) / (2 * dx)) * (phi(i + 1, j) - phi(i - 1, j)) -
					  (abs(u(i, j)) / (2 * dx)) * (phi(i + 1, j) - 2 * phi(i, j) + phi(i - 1, j)) +
					  (v(i, j) / (2 * dy)) * (phi(i, j + 1) - phi(i, j - 1)) -
					  (abs(v(i, j)) / 2 * dy) * (phi(i, j + 1) - 2 * phi(i, j) + phi(i, j - 1));
			*/
			//This one works because abs returns an int value 
			fc_Curr(i, j) = (-0.5 * sqrt(u(i, j) * u(i, j)) / dx) * (phi(i + 1, j) - 2 * phi(i, j) + phi(i - 1, j)) +
							(0.5 * u(i, j) / dx) * (phi(i + 1, j) - phi(i - 1, j)) -
							(0.5 * sqrt(v(i, j) * v(i, j)) / dy) * (phi(i, j + 1) - 2 * phi(i, j) + phi(i, j - 1)) +
							(0.5 * v(i, j) / dy) * (phi(i, j + 1) - phi(i, j - 1));

		}
}


void CDS2(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{

	unsigned long i, j;
	const double dx = G.dx();
	const double dy = G.dy();
	// Interior points only
	for(i = 1; i < G.Nx() - 1; i++)
		for(j = 1; j < G.Ny() - 1; j++)
		{
			fc_Curr(i, j) = (0.5 * u(i, j) / dx) * (phi(i + 1, j) - phi(i - 1, j)) +
								(0.5 * v(i, j) / dx) * (phi(i, j + 1) - phi(i, j - 1));
		}
}

void SOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
	unsigned long i, j;
	const double dx = G.dx();
	const double dy = G.dy();
	// Interior points only
	for(i = 1; i < G.Nx() - 1; i++)
		for(j = 1; j < G.Ny() - 1; j++)
		{	
			// Central Difference for points near boundary
			if(i == 1 || i == G.Nx() - 2 || j == 1 || j == G.Ny() - 2)
			{
				fc_Curr(i, j) = (0.5 * u(i, j) / dx) * (phi(i + 1, j) - phi(i - 1, j)) +
								(0.5 * v(i, j) / dx) * (phi(i, j + 1) - phi(i, j - 1));
			}
			// Second Order Upwind for Interior Points
			else
			{
				
				fc_Curr(i, j) = (0.25 * u(i, j) / dx) * (-4 * phi(i - 1, j) + phi(i - 2, j) + 4 * phi(i + 1, j) - phi(i + 2,j)) + 
					  (0.25 * sqrt(((u(i, j)) / dx) * ((u(i, j)) / dx))) * (6 * phi(i, j) - 4 * phi(i - 1, j) + phi(i - 2, j) - 4 * phi(i + 1, j) + phi(i + 2, j)) +
					  (0.25 * v(i, j) / dy) * (-4 * phi(i, j - 1) + phi(i, j - 2) + 4 * phi(i, j + 1) - phi(i, j + 2)) +  
					  (0.25 * sqrt(((v(i, j)) / dy) * ((v(i, j)) / dy))) * (6 * phi(i, j) - 4 * phi(i, j - 1) + phi(i, j - 2) - 4 * phi(i, j + 1) + phi(i, j + 2));
			}
		}
}

void eulerExp(Vector &phi, const Vector &fc_Curr, double &dt)
{
	phi = phi - (dt * fc_Curr);
}

void abs2Exp(Vector &phi, const Vector &fc_Curr, const Vector &fc_Prev, double &dt)
{
	phi = phi - dt * (1.5 * fc_Curr - 0.5 * fc_Prev);
}

// Add code from project 1
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
			Boundary conditions the same here, but with the difference (I/dt - 0.5 * A)
			I is 1 along the diagnol, and 0 everywhere else. Therefore if the current position is not along the diagnol I/dt becomes 0, 
			and if it is, it becomes 1/dt. Only M(i, j, 2) is along the diagnol.
			*/
			M(i,j,0) = -0.5 * (1 / (dx * dx));;
			M(i,j,1) = -0.5 * (1 / (dy * dy));
			M(i,j,2) = (1 / dt) - 0.5 * (- 2 / (dx * dx) - 2 / (dy * dy));
			M(i,j,3) = -0.5 * (1 / (dy * dy));
			M(i,j,4) = -0.5 * (1 / (dx * dx));
		}

	for(i = 0; i < Nx; i++)
		for(int t = 0; t < 5; t++){
			if(t == 2)
			{
				// Update A matrix for points we are updating
				M(i, 0, t) = (-1 / dy);
				M(i, Ny - 1, t) = (1 / dy);
			}
			else if(t == 3)
			{
				// Update A matrix for j + 1 points 
				M(i, 0, t) = (1 / dy);
				M(i, Ny - 1, t) = b;  // This point isnt used for the top boundary
			}
			else if(t == 1)
			{
				// Update A matrix for j - 1 points
				M(i, Ny - 1, t) = (-1 / dy);
				M(i, 0, t) = b;  // This point isnt used for bottom boundary 
			}
			else
			{
				// Set other points to be 0
				M(i, 0, t) = b;
				M(i, Ny - 1, t) = b; 
			}
		}

	for(j = 0; j < Ny; j++)
		for(int t = 0; t < 5; t++){
			M(0,j,t) = (t == 2 ? a : b);
			M(Nx-1,j,t) = (t == 2 ? a : b);
		}
}

void computeDiffusion(Vector &R, const Vector &phi, const Grid &G, double alpha)
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
			R(i,j) = alpha * (((phi(i-1, j) - (2 * phi(i, j)) + phi(i + 1, j)) / (dx * dx)) + 
					 ((phi(i, j-1) - (2 * phi(i, j)) + phi(i, j + 1)) / (dy * dy)));

}

void applyBC(Vector &R, Vector &dphi, const Vector &phi, const Vector &v, const Grid &G, double alpha)
{
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dy = G.dy();
	double h0; // Double to store bottom boundary condition
	double h1; // Double to store top boundary condition
	for(i = 0; i < Nx; i ++)
	{
		// Applying Neumann boundary condition to the bottom boundary
		j = 0;
		h0 = 0.1 * (-v(i, j) / alpha * (1 - exp(v(i, j) / alpha)));
		R(i, j) = h0 - ((phi(i, j + 1) - phi(i, j)) / dy);
		dphi(i, j) = 0.0;
		// Apply Neumann bounday condition to the top Boundary 
		j = (Ny - 1);
		h1 = 0.1 * (-v(i, j) * exp(v(i, j) / alpha) / (alpha * (1 - exp(v(i, j) / alpha))));
		R(i, j) = h1 - ((phi(i, j) - phi(i, j - 1)) / dy);
		dphi(i, j) = 0.0;
	}
	// Apply Dirichlet boundary condition to left and right
	for(j = 0; j < Ny; j++)
	{
		// Apply Dirichlet bboundary condition to the left boundary 
		i = 0;
		R(i, j) = dphi(i, j) = 0.0;
		// Apply Dirichlet boundary condition to the right boundary 
		i = Nx - 1;
		R(i, j) = dphi(i, j) = 0.0;
	}
}

void updateBoundaries(Vector &phi,const Vector &v, const Grid &G, double alpha)
{
	double dy = G.dy();
	unsigned long Ny = G.Ny();
	unsigned long Nx = G.Nx();
	unsigned long i, j;
	double h0;
	double h1;

	for(i = 1; i < Nx - 1; i++)
	{
		j = 0;
		h0 = 0.1 * (-v(i, j) / alpha * (1 - exp(v(i, j) / alpha))); // Prescribed value for Neumann condition on bottom boundary
		phi(i, j) = phi(i, j + 1) - h0 * dy;

		j = (Ny - 1);
		h1 = 0.1 * (-v(i, j) * exp(v(i, j) / alpha) / (alpha * (1 - exp(v(i, j) / alpha)))); // Prescribed value for Neumann condition on top boundary
		phi(i, j) = phi(i, j - 1) - h1 * dy;
	}
}



void SolveConvection(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme, double CFL)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	Vector phi(Nx, Ny);
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);
	InitializePhi(phi,G);
	InitializeVel(u,v,G);

	// fc_Curr: forward convection vector at current time step
	// fc_Prev: forward convection vector at previous time step

	Vector fc_Curr(Nx, Ny), fc_Prev(Nx,Ny);

	int itime = 0; double time = 0.0; bool last = false;
	while(time < tf)
	{
		// record the previous convection vector
		fc_Prev = fc_Curr;
		// calculate the convection vector at the current time step
		switch(conScheme){
			case 1:
				FOU(fc_Curr, phi, u, v, G);
				break;
			case 2:
				CDS2(fc_Curr, phi, u, v, G);
				break;
			case 3:
				SOU(fc_Curr, phi, u, v, G);
				break;
			default:
				printf("invalid convection scheme.\n");
				exit(0);
		}
		switch(timeScheme){
			case 1:
				eulerExp(phi, fc_Curr, dt);
				break;
			case 2:
				itime < 1 ? eulerExp(phi, fc_Curr, dt) : abs2Exp(phi, fc_Curr, fc_Prev, dt);
				break;
			default:
				printf("invalid time-integration scheme.\n");
				exit(0);
		}

		time = (last == true ? tf :(itime+1)*dt);
		if(itime % 1 == 0 || last)
			printf("Time = %lf, maxphi = %14.12e, minphi = %14.12e\n",time,phi.GetMax(),phi.GetMin());

		if(last)
			break;

		++itime;
		if(itime*dt > tf){
			dt = tf - (itime-1)*dt;
			last = true;
		}
	}

	// Save VTK file with timestep in it
	std::string str_cfl = std::to_string(CFL);
	std::string phi_dt = "phi_CFL=";
	std::string vtk = ".vtk";
	std::string total = phi_dt + str_cfl + vtk;
	char fname[total.length() + 1];
	strcpy(fname, total.c_str());
	//char fname[20] = "phi_dt=.vtk";
	storeVTKStructured(phi, G, fname);
}

void SolveConvectionDiffusion(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme, const double alpha)
{

	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	
	Vector phi(Nx, Ny);
	// Vector phi_conv(Nx, Ny); // Temporary vector to store convection phi values
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);
	Matrix A(Nx,Ny); // LHS matrix A
	// forward convection current and previous
	Vector fc_Curr(Nx,Ny);
	Vector fc_Prev(Nx,Ny);
	Vector R(Nx,Ny); // Residual R
	Vector id(Nx,Ny); // residual from implicite diffusion Aphi^n in the previous project
	Vector dphi(Nx,Ny); // Increment dphi
		

	//Initialize A, Phi, phi_conv and Vel
	computeTransientMatrix(A, G, dt);
	InitializeVelCD(u, v, G);
	InitializePhiCD(phi, u, v, G, alpha);
	// InitializePhiCD(phi_conv, u, v, G, alpha);

	// Compute diffusion and residual norm 
	computeDiffusion(id, phi, G, alpha);
	CDS2(fc_Curr, phi, u, v, G);
	R = id - 1.5 * fc_Curr + 0.5 * fc_Prev;  // Update residual to subtract convective solution
	applyBC(R, dphi, phi, v, G, alpha);
	double R0 = 0.0;
	R0 = R.L2Norm();
	printf("Initial Residual Norm %14.12e\n", R0);
	printf("Original dphi norm %14.12e\n", dphi.L2Norm());

	// Start time loop 
	int itime = 0;

    while(itime * dt <= tf)
	{
		// Record previous convection vector 
		fc_Prev = fc_Curr;
		solveGS(dphi, A, R);  // Solve Adphi = R
		phi = phi + dphi;  // Update phi
		updateBoundaries(phi, v, G, alpha);  // Re initialize phi points based on Neuman condition - they were overwritten in last line 
		computeDiffusion(id, phi, G, alpha);  // Solve for diffusion residual 
		// calculate the convection vector at the current time step
		switch(conScheme)
		{
			case 1:
				FOU(fc_Curr, phi, u, v, G);
				break;
			case 2:
				CDS2(fc_Curr, phi, u, v, G);
				break;
			case 3:
				SOU(fc_Curr, phi, u, v, G);
				break;
			default:
				printf("invalid convection scheme.\n");
				exit(0);
		}
		R = id - 1.5 * fc_Curr + 0.5 * fc_Prev;  // Update residual to include convection terms
		applyBC(R, dphi, phi, v, G, alpha); // Apply boundary conditions 
		double R1 = R.L2Norm();  // Calculate R1 to find residual norm 
		printf("Time-Step = %d\n",++itime); 
		printf("Convergence Criteria %12.12e\n Current Time = %lf\n", dphi.L2Norm(), itime*dt);
		//Check convergence
		if(dphi.L2Norm() < 1e-8)
		{
			printf("Steady state reached in %d time steps.\n Final time = %lf.\n",itime,itime*dt);
			break;
		}
	}
	char fname[20] = "Phi_2_working.vtk";
	storeVTKStructured(phi, G, fname);
}





int main()
{
	// method and problem
	unsigned short conScheme = 3;
	unsigned short timeScheme = 2;
	unsigned short pbtype = 1;

	switch (conScheme){
		case 1: std::cout << "First Order Upwind is used" << std::endl; break;
		case 2: std::cout << "Second order central difference is used" << std::endl; break;
		case 3: std::cout << "Second order upwind is used" << std::endl; break;
		default: std::cout << "invalid convection scheme" << std::endl; exit(0);
	}

	switch (timeScheme){
		case 1: std::cout << "Forward Euler is used" << std::endl; break;
		case 2: std::cout << "Adams-Bashforth 2 is used" << std::endl; break;
		default: std::cout << "invalid time-integration scheme" << std::endl; exit(0);
	}

	if (pbtype == 1)
	{
		std::cout << "2D convection equation is solved" << std::endl; 
		unsigned long Nx1, Ny1;
		double xlim1[2] = {-1, 1};
		double ylim1[2] = {-1, 1};
		Nx1 = Ny1 = 201;
		int numcheck = 9; // Number of CFL numbers to run
		double tf1 = 2*PI;
		double CFL[numcheck] = {0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195};
		// dt1 originally 0.05*2/200
		Grid G1(Nx1,Ny1,xlim1,ylim1);
		for(int i = 0; i < numcheck; i++)
		{
			double dt1 = CFL[i]*2.0/200.0;
			SolveConvection(G1, tf1, dt1, conScheme, timeScheme, CFL[i]);
		}
		

	} 
	else if (pbtype == 2)
	{
		std::cout << "2D convection-diffusion equation is solved" << std::endl;
		// Variables to define Grid
		unsigned long Nx2, Ny2;
		double xlim2[2] = {0, 1};
		double ylim2[2] = {0, 1};
		Nx2 = Ny2 = 17;
		// Define Problem
		double tf2 = 2 * 3.1415;
		//double tf2 = 1*0.05/200;
		double dt2 = 0.05/(Nx2 - 1);
		double alpha = 0.1;
		// Initialize and solve
		Grid G2(Nx2,Ny2,xlim2,ylim2);
		SolveConvectionDiffusion(G2, tf2, dt2, conScheme, timeScheme, alpha);
	}
	else if (pbtype == 3)
	{
		std::cout << "2D convection-diffusion exact equation is solved\n";
		// Initialize Grid
		unsigned long Nx3, Ny3;
	 	double xlim3[2] = {0, 1};
	 	double ylim3[2] = {0, 1};
	 	Nx3 = Ny3 = 65;
	 	Grid G3(Nx3, Ny3, xlim3, ylim3);
	 	// Define diffusivity constant
	 	double alpha = 0.1;

	 	// Solve exact equation
	 	InitializePhiCD_Exact(G3, alpha);

	}
	else 
	{
	 	std::cout << "Problem undefined" << std::endl;
	}
}