#include "Base.h"


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

void InitializePhiCD(Vector &phi, const Grid &G);
void InitializeVelCD(Vector &u, Vector &v, const Grid &G);

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

double f1(double x, double y){
	return 1.0 * exp(-1500 * (pow((x - 0.5), 2) + pow((y - 0.5), 2)));
}

void InitializePhiCD(Vector &phi, const Grid &G)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i,j;

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
			phi(i,j) = f1( G.x(i) , G.y(j) );
}

void InitializeVelCD(Vector &u, Vector &v, const Grid &G)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	unsigned long i, j;
	double x,y;

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++){
			u(i,j) = 1.0;
			v(i,j) = 0.05;
		}
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

void SolveConvection(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme)
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

	char fname[20] = "Phi.vtk";
	storeVTKStructured(phi, G, fname);


}

void SolveConvectionDiffusion(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme)
{

	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	
	Vector phi(Nx, Ny);
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);

	InitializePhiCD(phi,G);
	InitializeVelCD(u,v,G);

	// variables needed for Newton-Raphson iteration
	// LHS matrix A
	Matrix A(Nx,Ny);
	// Residual R
	Vector R(Nx,Ny);
	// Increment dphi
	Vector dphi(Nx,Ny);
	
	// variables needed for calculating the residual
	// forward convection current and previous
	Vector fc_Curr(Nx,Ny);
	Vector fc_Prev(Nx,Ny);
	// residual from implicite diffusion Aphi^n in the previous project
	Vector id(Nx,Ny);

	
}





int main()
{
	// method and problem
	unsigned short conScheme = 2;
	unsigned short timeScheme = 1;
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

	if (pbtype == 1){

		std::cout << "2D convection equation is solved" << std::endl; 
		unsigned long Nx1, Ny1;
		double xlim1[2] = {-1, 1};
		double ylim1[2] = {-1, 1};
		Nx1 = Ny1 = 201;
		double tf1 = 2*PI;
		double dt1 = 0.2*2/200;
		// dt1 originally 0.05*2/200
		Grid G1(Nx1,Ny1,xlim1,ylim1);
		SolveConvection(G1, tf1, dt1, conScheme, timeScheme);

	} else if (pbtype == 2){

		std::cout << "2D convection-diffusion equation is solved" << std::endl;
		unsigned long Nx2, Ny2;
		double xlim2[2] = {0, 1};
		double ylim2[2] = {0, 1};
		Nx2 = Ny2 = 101;
		double tf2 = 2*PI;
		double dt2 = 0.05*2/200;
		Grid G2(Nx2,Ny2,xlim2,ylim2);
		SolveConvectionDiffusion(G2, tf2, dt2, conScheme, timeScheme);

	} else {
	 	std::cout << "Problem undefined" << std::endl;
	}
}
