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
void ResetPhi(Vector& phi, const Grid& G);
void InitializeVel(Vector &u, Vector &v, const Grid &G);
void ExactPhi(Vector& phi, Vector&u, Vector&v, const Grid&G);

/*================================================================================================
 * Correct the initialization function for the 2D convection diffusion equation f1 
 * pay attention to boundary values (use the value from the boundary condition)
 *================================================================================================
 */
double f1(double x, double y);
/*================================================================================================*/

void InitializePhiCD(Vector &phi, const Vector& u, const Vector& v, const Grid &G, const double &Re);
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
void SolveConvectionDiffusionX(const Grid &G, double tf, double dt, const unsigned short scheme, double Re);
void SolveConvectionDiffusionY(const Grid& G, double tf, double dt, const unsigned short scheme, double Re);
/*================================================================================================*/

void computeTransientMatrix(Matrix& M, const Grid& G, const double& dt, const double& Re);
void computeDiffusion(Vector& R, const Vector& phi, const Vector& u, const Vector& v, const Grid& G);
void applyBC(Vector& R, Vector& dphi, const Vector& phi, const Vector& u, const Vector& v, const Grid& G, const double &Re);

double f0(double x, double y){
	return 5.0 * exp(-1500.0 * ((x - 0.25) * (x - 0.25) + y * y));
}

double f1(double x, double y) {
	return 1.0 * exp(-1500 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
}
	
void InitializePhi(Vector& phi, const Grid& G)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i, j;
	double x, y;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			phi(i, j) = f0(G.x(i), G.y(j));
}

void ResetPhi(Vector& phi, const Grid& G)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i, j;
	double x, y;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			phi(i, j) = 0;

}

void ExactPhi2(Vector& phi_exact, const Vector& u, const Vector& v, const Grid& G, const double& Re)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i, j;
	double dx = G.dx();
	double dy = G.dy();
	double a = Re;

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			phi_exact(i, j) = 5.0 * ((1 - exp(u(i, j) * i * dx / a)) / (1 - exp(u(i, j) / a))) + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));
		}
	}
	char fname[20] = "Phi_Exact2.vtk";
	storeVTKStructured(phi_exact, G, fname);
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

void InitializePhiCD(Vector &phi, const Vector& u, const Vector& v, const Grid &G, const double& Re)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i,j;
	double dx = G.dx();
	double dy = G.dy();
	const double a = 1 / Re;
	long long phi_long;

	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Ny; j++) {
			phi(i, j) = f1(G.x(i), G.y(j));
		}
	}
	/*
	for (i = 0; i < Nx; i++) {
		j = 0;
		phi(i, j) = 0.1 * (-v(i, j) / (a * (1 - exp(v(i, j) / a))));    //This is the bottom boundary

		j = Ny - 1;
		phi(i, j) = 0.1 * (-v(i, j) * exp(v(i, j) / a) / (a * (1 - exp(v(i, j) / a))));    //This is the top boundary
	}
	*/
	for (j = 0; j < Ny; j++) {
		i = 0;
		phi(i, j) = 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));    // This is the left boundary

		i = Nx - 1;
		phi(i, j) = 5.0 + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));        //This is the right boundary

	}

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
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();
	double u_plus, u_minus, v_plus, v_minus, u_1, u_2, v_1, v_2;

	for (i = 1; i < Nx - 1; i++) {
		for (j = 1; j < Ny - 1; j++) {

			u_plus = 0.5 * (u(i, j) + pow(pow(u(i, j), 2), 0.5));
			u_minus = 0.5 * (u(i, j) - pow(pow(u(i, j), 2), 0.5));
			v_plus = 0.5 * (v(i, j) + pow(pow(v(i, j), 2), 0.5));
			v_minus = 0.5 * (v(i, j) - pow(pow(v(i, j), 2), 0.5));
			
			u_1 = (phi(i, j) - phi(i - 1, j)) / dx;
			u_2 = (phi(i + 1, j) - phi(i, j)) / dx;
			v_1 = (phi(i, j) - phi(i, j - 1)) / dy;
			v_2 = (phi(i, j + 1) - phi(i, j)) / dy;
			fc_Curr(i, j) = -u_plus * u_1 - u_minus * u_2 - v_plus * v_1 - v_minus * v_2;
			//std::cout << u_plus << u_minus << v_plus << v_minus << std::endl;
		}
	}
}

void CDS2(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
	unsigned long i, j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();
	double u_1, v_1;

	for (i = 1; i < Nx - 1; i++) {
		for (j = 1; j < Ny - 1; j++) {
			u_1 = u(i, j) * ((phi(i + 1, j) - phi(i - 1, j)) / (2.0 * dx));
			v_1 = v(i, j) * ((phi(i, j + 1) - phi(i, j - 1)) / (2.0 * dy));
			fc_Curr(i, j) = -u_1 - v_1;
			//std::cout << u_plus << u_minus << v_plus << v_minus << std::endl;
		}
	}
}

void SOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
	unsigned long i, j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();
	double u_plus, u_minus, v_plus, v_minus, u_1, u_2, v_1, v_2, u_cd, v_cd;

	for (i = 1; i < Nx - 1; i++) {
		for (j = 1; j < Ny - 1; j++) {
			if (i == 1 || i == Nx - 1 || j == 1 || j == Ny - 1) {
				// Do Central Difference for the nodes next to the boundary
				u_cd = u(i, j) * ((phi(i + 1, j) - phi(i - 1, j)) / (2.0 * dx));
				v_cd = v(i, j) * ((phi(i, j + 1) - phi(i, j - 1)) / (2.0 * dy));
				fc_Curr(i, j) = -u_cd - v_cd;
			}

			else {
				// Do SOU for the rest of the nodes
				u_plus = 0.5 * (u(i, j) + pow(pow(u(i, j), 2), 0.5));
				u_minus = 0.5 * (u(i, j) - pow(pow(u(i, j), 2), 0.5));
				v_plus = 0.5 * (v(i, j) + pow(pow(v(i, j), 2), 0.5));
				v_minus = 0.5 * (v(i, j) - pow(pow(v(i, j), 2), 0.5));
				
				u_1 = (3 * phi(i, j) - 4 * phi(i - 1, j) + phi(i - 2, j)) / (2.0 * dx);
				u_2 = (-3 * phi(i, j) + 4 * phi(i + 1, j) - phi(i + 2, j)) / (2.0 * dx);
				v_1 = (3 * phi(i, j) - 4 * phi(i, j - 1) + phi(i, j - 2)) / (2.0 * dy);
				v_2 = (-3 * phi(i, j) + 4 * phi(i, j + 1) - phi(i, j + 2)) / (2.0 * dy);
				fc_Curr(i, j) = -u_plus * u_1 - u_minus * u_2 - v_plus * v_1 - v_minus * v_2;

				//std::cout << u_plus << u_minus << v_plus << v_minus << std::endl;
			}
		}
	}

}

void eulerExp(Vector &phi, const Vector &fc_Curr, double &dt)
{
	phi = phi + dt * fc_Curr;
}

void abs2Exp(Vector &phi, const Vector &fc_Curr, const Vector &fc_Prev, double &dt)
{
	phi = phi + dt * (1.5 * fc_Curr - 0.5 * fc_Prev);
}

void computeTransientMatrix(Matrix& M, const Grid& G, const double& dt, const double& Re)
{
	unsigned long i, j;
	const double dx = G.dx();
	const double dy = G.dy();
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();
	
	// a Trapezoidal
	const double a = 1.0;
	const double b = 0;

	for (i = 1; i < Nx - 1; i++) {
		for (j = 1; j < Ny - 1; j++) {
			// Trapezoidal Rule for temporal integration
			M(i, j, 0) = -1.0 / 2.0 / Re * (1.0 / dx / dx);
			M(i, j, 1) = -1.0 / 2.0 / Re * (1.0 / dy / dy);
			M(i, j, 3) = -1.0 / 2.0 / Re * (1.0 / dy / dy);
			M(i, j, 4) = -1.0 / 2.0 / Re * (1.0 / dx / dx);
			M(i, j, 2) = 1.0 / dt - 1.0 / 2.0 / Re * (-2.0 / dx / dx - 2.0 / dy / dy);
			//M(i, j, 2) = -(M(i, j, 0) + M(i, j, 1) + M(i, j, 3) + M(i, j, 4)) + 1.0 / dt;
		}
	}

	// Need to update with Neumann BC in Matrix (-1 values)
	// These are outer perimeter values for A matrix
	for (i = 0; i < Nx; i++) {
		for (int t = 0; t < 5; t++) {
			M(i, 0, t) = (t == 2 ? -1.0 : b);     // This is bottom boundary
			M(i, Ny - 1, t) = (t == 2 ? 1.0 : b);     // This is top boundary
		}
	}

	for (j = 0; j < Ny; j++) {
		for (int t = 0; t < 5; t++) {
			M(0, j, t) = (t == 2 ? a : b);     // This is left boundary
			M(Nx - 1, j, t) = (t == 2 ? a : b);     // This is right boundary
		}
	}
}

void computeDiffusion(Vector& R, const Vector& phi, const Vector& u, const Vector& v, const Grid& G)
{
	unsigned long i, j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();
	//const double a = 0.1; //Make this the diffusion coefficient alpha

	// This is RHS equation
	// Add the other 2 Neumann BC conditions here
	for (i = 1; i < Nx - 1; i++) {
		for (j = 1; j < Ny - 1; j++) {
			R(i, j) = (phi(i - 1, j) - 2 * phi(i, j) + phi(i + 1, j)) / (dx * dx) + (phi(i, j - 1) - 2 * phi(i, j) + phi(i, j + 1)) / (dy * dy);
		}
	}
	/*
	for (i = 0; i < Nx; i++) {
		j = 0;
		R(i, j) = (dy * (0.1 * (-v(i, j) / (a * (1 - exp(v(i, j) / a))))) - phi(i, 0) + phi(i, 1));     //This is the bottom boundary

		j = Ny - 1;
		R(i, j) = (dy * (0.1 * (-v(i, j) * exp(v(i, j) / a) / (a * (1 - exp(v(i, j) / a))))) - phi(i, j) + phi(i, j - 1));     //This is the top boundary
	}
	*/
}

void applyBC(Vector& R, Vector& dphi, const Vector& phi, const Vector& u, const Vector& v, const Grid& G, const double& Re)
{
	unsigned long i, j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dy = G.dy();
	double dx = G.dx();

	// This does outer boundary conditions

	const double a = 1 / Re;

	// Distance in y-direction = j*dy

	for (i = 0; i < Nx; i++) {
		j = 0;
		R(i, j) = dphi(i, j) = (dy * (0.1 * (-v(i, j) / (a * (1 - exp(v(i, j) / a))))) + phi(i, 0) - phi(i, 1));     //This is the bottom boundary
		//R(i, j) = (dy * (0.1 * (-v(i, j) / (a * (1 - exp(v(i, j) / a))))) + phi(i, 0) - phi(i, 1));     //This is the bottom boundary
		//R(i, j) = dphi(i, j) = 1 / dy * (dy * (0.1 * (-v(i, j) / (a * (1 - exp(v(i, j) / a))))) + phi(i, 0) - phi(i, 1));
		//R(i, j) = 0.1 * (-v(i, j) / (a * (1 - exp(v(i, j) / a))));
		//R(i, j) = dphi(i, j) = 0.1 * (-v(i, j) / (a * (1 - exp(v(i, j) / a))));    //This is the bottom boundary
		//R(i, j) = dphi(i, j) = 1.0 * exp(-1500.0 * ((i * dx - 0.5) * (i * dx - 0.5) + (j * dy - 0.5) * (j * dy - 0.5)));    // This is the exact solution for the bottom boundary
		//R(i, j) = dphi(i, j) = 0;
		//R(i, j) = dphi(i, j) = 5.0 * ((1 - exp(u(i, j) * i * dx / a)) / (1 - exp(u(i, j) / a))) + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));
		//R(i, j) = 5.0 * ((1 - exp(u(i, j) * i * dx / a)) / (1 - exp(u(i, j) / a))) + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));
		//dphi(i, j) = (dy * (0.1 * (-v(i, j) / (a * (1 - exp(v(i, j) / a))))) + phi(i, 0) - phi(i, 1));

		j = Ny - 1;
		R(i, j) = dphi(i, j) = (dy * (0.1 * (-v(i, j) * exp(v(i, j) / a) / (a * (1 - exp(v(i, j) / a))))) - phi(i, j) + phi(i, j - 1));     //This is the top boundary
		//R(i, j) = (dy * (0.1 * (-v(i, j) * exp(v(i, j) / a) / (a * (1 - exp(v(i, j) / a))))) - phi(i, j) + phi(i, j - 1));     //This is the top boundary
		//dphi(i, j) = 10;
		//R(i, j) = dphi(i, j) = 1 / dy * (dy * (0.1 * (-v(i, j) * exp(v(i, j) / a) / (a * (1 - exp(v(i, j) / a))))) - phi(i, j) + phi(i, j - 1));    //This is the top boundary
		//R(i, j) = 0.1 * (-v(i, j) * exp(v(i, j) / a) / (a * (1 - exp(v(i, j) / a))));     //This is the top boundary
		//R(i, j) = dphi(i, j) = 0.1 * (-v(i, j) * exp(v(i, j) / a) / (a * (1 - exp(v(i, j) / a))));    //This is the top boundary
		//R(i, j) = dphi(i, j) = 1.0 * exp(-1500.0 * ((i * dx - 0.5) * (i * dx - 0.5) + (j * dy - 0.5) * (j * dy - 0.5)));    // This is the exact solution for the top boundary
		//R(i, j) = dphi(i, j) = 0;
		//R(i, j) = dphi(i, j) = 5.0 * ((1 - exp(u(i, j) * i * dx / a)) / (1 - exp(u(i, j) / a))) + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));
		//R(i, j) = 5.0 * ((1 - exp(u(i, j) * i * dx / a)) / (1 - exp(u(i, j) / a))) + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));
		//dphi(i, j) = (dy * (0.1 * (-v(i, j) * exp(v(i, j) / a) / (a * (1 - exp(v(i, j) / a))))) - phi(i, j) + phi(i, j - 1));
	}

	for (j = 0; j < Ny; j++) {
		i = 0;
		//R(i, j) = dphi(i, j) = 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));    // This is the left boundary
		//R(i, j) = 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));
		//dphi(i, j) = 0;
		//R(i, j) = dphi(i, j) = 1.0 * exp(-1500.0 * ((i * dx - 0.5) * (i * dx - 0.5) + (j * dy - 0.5) * (j * dy - 0.5)));    // This is the exact solution for the left boundary
		R(i, j) = dphi(i, j) = 0;     // This is the left boundary
		//R(i, j) = dphi(i, j) = 5.0 * ((1 - exp(u(i, j) * i * dx / a)) / (1 - exp(u(i, j) / a))) + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));

		i = Nx - 1;
		//R(i, j) = dphi(i, j) = 5.0 + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));        //This is the right boundary
		//R(i, j) = 5.0 + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)))	;
		//dphi(i, j) = 0;
		//R(i, j) = 5.0 + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));        //This is the right boundary
		//R(i, j) = dphi(i, j) = 1.0 * exp(-1500.0 * ((i * dx - 0.5) * (i * dx - 0.5) + (j * dy - 0.5) * (j * dy - 0.5)));    // This is the exact solution for the right boundary
		R(i, j) = dphi(i, j) = 0;     // This is the right boundary
		//R(i, j) = dphi(i, j) = 5.0 * ((1 - exp(u(i, j) * i * dx / a)) / (1 - exp(u(i, j) / a))) + 0.1 * ((1 - exp(v(i, j) * j * dy / a)) / (1 - exp(v(i, j) / a)));
	}
}

void SolveConvection(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	Vector phi(Nx,Ny);
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);
	InitializePhi(phi, G);
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

void SolveConvectionDiffusionX(const Grid& G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme, double& Re)
{
	const double exact = 1;
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	Vector phi(Nx, Ny);
	Vector phi_exact(Nx, Ny);
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);

	Matrix A(Nx, Ny);

	// Residual Vector
	Vector R(Nx, Ny);
	Vector R_diff(Nx, Ny);

	// Increment dphi
	Vector dphi(Nx, Ny);
	Vector id(Nx, Ny);

	InitializeVelCD(u, v, G);
	InitializePhiCD(phi, u, v, G, Re);

	ExactPhi2(phi_exact, u, v, G, Re);

	// fc_Curr: forward convection vector at current time step
	// fc_Prev: forward convection vector at previous time step

	Vector fc_Curr(Nx, Ny), fc_Prev(Nx, Ny);

	//Compute Transient (LHS) Matrix = (I/dt - 0.5*D) from Project 1 with new BC's
	computeTransientMatrix(A, G, dt, Re);

	//Compute Residual and Residual norm with new BC's
	double R0 = 0.0;
	computeDiffusion(R_diff, phi, u, v, G);
	CDS2(fc_Curr, phi, u, v, G);

	//R = R_diff - 1.5 * fc_Curr + 0.5 * fc_Prev;
	R = (R_diff + 1.5 * fc_Curr - 0.5 * fc_Prev);
	//R = R_diff;

	applyBC(R, dphi, phi, u, v, G, Re);
	R0 = R.L2Norm();
	//printf("Initial Residual Norm = %14.12e\n", R0);

	int itime = 0; double time = 0.0; bool last = false;
	while (time < tf) {

		fc_Prev = fc_Curr;
		//ResetPhi(phi, G);

		//Solve the linear system Adu = R using Gauss Seidal (NR loop)
		solveGS(dphi, A, R);

		//Update the solution phi
		phi = phi + dphi;

		// Update R
		computeDiffusion(R_diff, phi, u, v, G);
		CDS2(fc_Curr, phi, u, v, G);
		//R = R_diff - 1.5 * fc_Curr + 0.5 * fc_Prev;
		R = (R_diff + 1.5 * fc_Curr - 0.5 * fc_Prev);
		//R = R_diff;

		// Solve new BC's
		applyBC(R, dphi, phi, u, v, G, Re);
		double R1 = R.L2Norm();

		//Print timestep and residual norm
		//printf("Time-Step = %d\n", ++itime);
		//printf("Residual Norm = %14.12e,\n Residual Norm Ratio (R/R0) = %14.12e\n", R1, R1 / R0);
		if (itime % 100 == 1)
		{
			printf("Time-Step = %d\n", ++itime);
			printf("Residual Norm = %14.12e,\n Residual Norm Ratio (R/R0) = %14.12e\n", R1, R1 / R0);
			printf("Time = %lf, maxphi = %14.12e, minphi = %14.12e\n", time, phi.GetMax(), phi.GetMin());
		}

		//Check convergence
		if (dphi.L2Norm() < 1e-8) {
			printf("Steady state reached in %d time steps.\n Final time = %lf.\n", itime, itime * dt);
			break;
		}

		//Print max and min phi
		time = (last == true ? tf : (itime + 1) * dt);
		if (itime % 1 == 0 || last)
			//printf("Time = %lf, maxphi = %14.12e, minphi = %14.12e\n", time, phi.GetMax(), phi.GetMin());

		if (last)
			break;

		//Update time
		++itime;
		if (itime * dt > tf) {
			dt = tf - (itime - 1) * dt;
			last = true;
		}

	}

	//Store solution in v_x vector
	Vector v_x(Nx, Ny);
	v_x = phi;

	char fname[20] = "v_x.vtk";
	storeVTKStructured(v_x, G, fname);

	//Compute difference error
	Vector e(Nx, Ny);
	unsigned long i, j;
	e = v_x - phi_exact;
	double Linf = e.LinfNorm(i, j);
	printf("Difference Error (L2, L_exact) = %14.12e, %14.12e at (%lu, %lu)\n", e.L2Norm(), Linf, i, j);
}


void SolveConvectionDiffusionY(const Grid& G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme, double Re)
{	
	const double exact = 1;
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	Vector phi(Nx, Ny);
	Vector phi_exact(Nx, Ny);
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);

	Matrix A(Nx, Ny);

	// Residual Vector
	Vector R(Nx, Ny);
	Vector R_diff(Nx, Ny);

	// Increment dphi
	Vector dphi(Nx, Ny);
	Vector id(Nx, Ny);

	InitializeVelCD(u, v, G);
	InitializePhiCD(phi, u, v, G, Re);

	ExactPhi2(phi_exact, u, v, G, Re);

	// fc_Curr: forward convection vector at current time step
	// fc_Prev: forward convection vector at previous time step

	Vector fc_Curr(Nx, Ny), fc_Prev(Nx, Ny);

	//Compute Transient (LHS) Matrix = (I/dt - 0.5*D) from Project 1 with new BC's
	computeTransientMatrix(A, G, dt, Re);

	//Compute Residual and Residual norm with new BC's
	double R0 = 0.0;
	computeDiffusion(R_diff, phi, u, v, G);
	CDS2(fc_Curr, phi, u, v, G);

	//R = R_diff - 1.5 * fc_Curr + 0.5 * fc_Prev;
	R = (R_diff + 1.5 * fc_Curr - 0.5 * fc_Prev);
	//R = R_diff;

	applyBC(R, dphi, phi, u, v, G, Re);
	R0 = R.L2Norm();
	// printf("Initial Residual Norm = %14.12e\n", R0);

	int itime = 0; double time = 0.0; bool last = false;
	while (time < tf) {

		fc_Prev = fc_Curr;
		//ResetPhi(phi, G);

		//Solve the linear system Adu = R using Gauss Seidal (NR loop)
		solveGS(dphi, A, R);

		//Update the solution phi
		phi = phi + dphi;

		// Update R
		computeDiffusion(R_diff, phi, u, v, G);
		CDS2(fc_Curr, phi, u, v, G);
		//R = R_diff - 1.5 * fc_Curr + 0.5 * fc_Prev;
		R = (R_diff + 1.5 * fc_Curr - 0.5 * fc_Prev);
		//R = R_diff;

		// Solve new BC's
		applyBC(R, dphi, phi, u, v, G, Re);
		double R1 = R.L2Norm();

		//Print timestep and residual norm
		//printf("Time-Step = %d\n", ++itime);
		//printf("Residual Norm = %14.12e,\n Residual Norm Ratio (R/R0) = %14.12e\n", R1, R1 / R0);
		if (itime % 100 == 1)
		{
			printf("Time-Step = %d\n", ++itime);
			printf("Residual Norm = %14.12e,\n Residual Norm Ratio (R/R0) = %14.12e\n", R1, R1 / R0);
			printf("Time = %lf, maxphi = %14.12e, minphi = %14.12e\n", time, phi.GetMax(), phi.GetMin());
		}

		//Check convergence
		if (dphi.L2Norm() < 1e-8) {
			printf("Steady state reached in %d time steps.\n Final time = %lf.\n", itime, itime * dt);
			break;
		}

		//Print max and min phi
		time = (last == true ? tf : (itime + 1) * dt);
		if (itime % 1 == 0 || last)
			//printf("Time = %lf, maxphi = %14.12e, minphi = %14.12e\n", time, phi.GetMax(), phi.GetMin());

		if (last)
			break;

		//Update time
		++itime;
		if (itime * dt > tf) {
			dt = tf - (itime - 1) * dt;
			last = true;
		}

	}

	//Store solution in v_y vector
	Vector v_y(Nx, Ny);
	v_y = phi;

	char fname[20] = "v_y.vtk";
	storeVTKStructured(v_y, G, fname);

	//Compute difference error
	Vector e(Nx, Ny);
	unsigned long i, j;
	e = v_y - phi_exact;
	double Linf = e.LinfNorm(i, j);
	printf("Difference Error (L2, L_exact) = %14.12e, %14.12e at (%lu, %lu)\n", e.L2Norm(), Linf, i, j);
}


int main()
{
	// method and problem
	unsigned short conScheme = 2;
	unsigned short timeScheme = 2;
	unsigned short pbtype = 2;

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
		double tf1 = 2*3.1416;
		//double dt1 = 0.05*2/200;
		//double dt1 = 2*0.2/200;
		double dt1 = 0.28 * 2.0 / (Nx1 - 1);
		Grid G1(Nx1,Ny1,xlim1,ylim1);
		SolveConvection(G1, tf1, dt1, conScheme, timeScheme);

	} else if (pbtype == 2){
		double Re;
		Re = 100;
	
		std::cout << "2D convection-diffusion equation is solved" << std::endl;
		unsigned long Nx2, Ny2;
		double xlim2[2] = {0, 1};
		double ylim2[2] = {0, 1};
		Nx2 = Ny2 = 33;
		double tf2 = 2*3.1416;
		//double dt2 = 0.05 * 2 / Nx2;
		double dt2 = 0.99999999 / 500 / (Nx2 - 1);
		//double dt2 = 0.001;
		Grid G2(Nx2,Ny2,xlim2,ylim2);
		printf("======================================Solving Convection Diffusion in X-Direction======================================\n");
		SolveConvectionDiffusionX(G2, tf2, dt2, conScheme, timeScheme, Re);
		printf("======================================Solving Convection Diffusion in Y-Direction======================================\n");
		SolveConvectionDiffusionY(G2, tf2, dt2, conScheme, timeScheme, Re);

	} else {
	 	std::cout << "Problem undefined" << std::endl;
	}



}