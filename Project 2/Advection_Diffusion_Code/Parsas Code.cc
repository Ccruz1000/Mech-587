//applyBC
//applyDiffusivityalpha

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

void InitializePhiCD(Vector &u, Vector &v, Vector &phi, const Grid &G);
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
	return 5.0*exp(-1500.0*(pow((x-0.25),2)+pow(y,2)));
	//return 5.0*exp(-1500.0*((x-0.25)*(x-0.25)+y*y));
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

	for(i = 0; i < Nx; i++){
		for(j = 0; j < Ny; j++){
			u(i,j) = G.y(j);
			v(i,j) = -G.x(i);
		}
	}
}

double f1(double x, double y){
	return 1.0*exp(-1500.0*(pow((x-0.5),2)+pow((y-0.5),2)));
}
/*
void InitializePhiCD(Vector &u, Vector &v, Vector &phi, const Grid &G, double alpha)
{
	const double dx=G.dx();
	const double dy=G.dy();
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i,j;

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{			
			if (i == 0)
			{
				phi(i,j) = 0.1*(1.0-exp(v(i,j)*G.y(j)/alpha))/(1.0-exp(v(i,j)/alpha));
			}
			else if (i == Nx-1)
			{
				phi(i,j) = 5.0 + 0.1*(1.0-exp(v(i,j)*G.y(j)/alpha))/(1.0-exp(v(i,j)/alpha));
			}
			else if (j == 0)
			{
				phi(i,j) = (5.0*(1.0-exp(u(i,j)*G.x(i)/alpha))/(1.0-exp(u(i,j)/alpha)) + 
					        0.1*(1.0-exp(v(i,j)*G.y(j)/alpha))/(1.0-exp(v(i,j)/alpha)));
			}
			else if (j == Ny-1)
			{
				phi(i,j) = (5.0*(1.0-exp(u(i,j)*G.x(i)/alpha))/(1.0-exp(u(i,j)/alpha)) + 
						    0.1*(1.0-exp(v(i,j)*G.y(j)/alpha))/(1.0-exp(v(i,j)/alpha)));
			}
			else
			{
				phi(i,j) = f1( G.x(i) , G.y(j) );
			}
		}		
}
*/
void InitializePhiCD(const Vector &u, const Vector &v, Vector &phi, const Grid &G, double alpha)
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
			else if(j == 0 || j == (Ny - 1))
			{
				phi(i, j) = 5.0 * ((1 - exp(G.x(i) * u(i, j) / alpha)) / (1 - exp(u(i, j) / alpha))) + 
				            0.1 * ((1 - exp(G.y(j) * v(i, j) / alpha)) / (1 - exp(v(i, j) / alpha)));
		
			}
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

	for(i = 0; i < Nx; i++){
		for(j = 0; j < Ny; j++){
			u(i,j) = 1.0;
			v(i,j) = 0.05;
		}
	}
}

void FOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
	const double dx=G.dx();
	const double dy=G.dy();
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i,j;

	for(i = 1; i < Nx-1; i++){
		for(j = 1; j < Ny-1; j++){			
			fc_Curr(i,j) = -0.5*(u(i,j))/dx*(phi(i+1,j)-phi(i-1,j)) 
					       +0.5*sqrt(((u(i,j))/dx)*((u(i,j))/dx))*(phi(i+1,j)-2*phi(i,j)+phi(i-1,j))
					       -0.5*(v(i,j))/dy*(phi(i,j+1)-phi(i,j-1)) 
					       +0.5*sqrt(((v(i,j))/dy)*((v(i,j))/dy))*(phi(i,j+1)-2*phi(i,j)+phi(i,j-1));
		}
	}
}


void CDS2(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
	const double dx=G.dx();
	const double dy=G.dy();
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i,j;

	for(i = 1; i < Nx-1; i++){
		for(j = 1; j < Ny-1; j++){	
			fc_Curr(i,j) = -0.5*u(i,j)/dx*(phi(i+1,j)+phi(i-1,j)-2*phi(i,j))
					       -0.5*v(i,j)/dy*(phi(i,j+1)+phi(i,j-1)-2*phi(i,j));
		}
	}
}

void SOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
	const double dx=G.dx();
	const double dy=G.dy();
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	unsigned long i,j;	
	for(i = 1; i < Nx-1; i++){
		for(j = 1; j < Ny-1; j++){	
			if (i == 1|| i == Nx-2 || j == 1 || j == Ny-2)
			{
			fc_Curr(i,j) = -0.5*u(i,j)/dx*(phi(i+1,j)+phi(i-1,j)-2*phi(i,j))
					       -0.5*v(i,j)/dy*(phi(i,j+1)+phi(i,j-1)-2*phi(i,j));
			}
					
			fc_Curr(i,j) = -0.25*u(i,j)/dx*(-4*phi(i-1,j)+phi(i-2,j)+4*phi(i+1,j)-phi(i+2,j)) 
					       -0.25*sqrt(((u(i,j))/dx)*((u(i,j))/dx))*(6*phi(i,j)-4*phi(i-1,j)+phi(i-2,j)-4*phi(i+1,j)+phi(i+2,j))
					       -0.25*v(i,j)/dy*(-4*phi(i,j-1)+phi(i,j-2)+4*phi(i,j+1)-phi(i,j+2)) 
					       -0.25*sqrt(((v(i,j))/dy)*((v(i,j))/dy))*(6*phi(i,j)-4*phi(i,j-1)+phi(i,j-2)-4*phi(i,j+1)+phi(i,j+2));
		}
	}

}

void eulerExp(Vector &phi, const Vector &fc_Curr, double &dt)
{
	phi=phi+fc_Curr*dt;
}

void abs2Exp(Vector &phi, const Vector &fc_Curr, const Vector &fc_Prev, double &dt)
{
	phi=(phi+dt/2*(3*fc_Curr-fc_Prev));
}

void SolveConvection(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	Vector phi(Nx,Ny);
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

/*============================================================================
 * This is the function for constructing the LHS matrix for transient solution

 *============================================================================
*/

void computeTransientMatrix(Matrix &M, const Grid &G, const double &dt)
{
	unsigned long i,j;
	const double dx = G.dx();
	const double dy = G.dy();
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	const double a = 1.0;
	const double b = 0.0;
	
	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++) {

			M(i,j,0) = -0.5*(1.0/(dx*dx));
			M(i,j,1) = -0.5*(1.0/(dy*dy));
			M(i,j,2) = 1.0/dt-0.5*(-2.0/(dx*dx)-2.0/(dy*dy));
			M(i,j,3) = -0.5*(1.0/(dy*dy));
			M(i,j,4) = -0.5*(1.0/(dx*dx));
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

 *===================================================================================
*/

void computeDiffusion(Vector &R, const Vector &phi, const Grid &G, double alpha)
{
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();


	for(i = 1; i < Nx-1; i++){
		for(j = 1; j < Ny-1; j++)
		{
			R(i,j) = alpha*((phi(i-1,j)-2.0*phi(i,j)+phi(i+1,j)) / (dx*dx) + 
				(phi(i,j-1)-2.0*phi(i,j)+phi(i,j+1)) / (dy*dy));

		}
	}

}


/*===================================================================================
 * This is the function for imposing bounadry condition
 * The variables you need to correct is:
 * a
 *===================================================================================
*/

void applyBC(Vector &R, Vector &dphi, const Grid &G)
{
 	const double dx=G.dx();
 	const double dy=G.dy();
 	unsigned long i,j;
 	unsigned long Nx = G.Nx();
 	unsigned long Ny = G.Ny();


 	for(i = 0; i < Nx; i++)
 	{
 		j = 0;
 		R(i,j) = dphi(i,j) = 0.0;//(5.0*(1.0-exp(u(i,j)*G.x(i)/alpha))/(1.0-exp(u(i,j)/alpha)) + 0.1*(1.0-exp(v(i,j)*G.y(j)/alpha))/(1.0-exp(v(i,j)/alpha)));//dy*0.1*-v(i,j)/(alpha*(1-exp(v(i,j)/alpha)))+dphi(i,j);//
 		j = (Ny-1);
 		R(i,j) = dphi(i,j) = 0.0;//(5.0*(1.0-exp(u(i,j)*G.x(i)/alpha))/(1.0-exp(u(i,j)/alpha)) + 0.1*(1.0-exp(v(i,j)*G.y(j)/alpha))/(1.0-exp(v(i,j)/alpha)));//dy*0.1*-v(i,j)*exp(v(i,j)/alpha)/(alpha*(1-exp(v(i,j)/alpha)))+dphi(i,j-1);//

 	}

 	for(j = 0; j < Ny; j++){
 		i = 0;
 		R(i,j) = dphi(i,j) = 0.0;//(0.1*(1.0-exp(v(i,j)*G.y(j)/alpha))/(1.0-exp(v(i,j)/alpha)));
 		i = (Nx-1);
 		R(i,j) = dphi(i,j) = 0.0;//(5.0 + 0.1*(1.0-exp(v(i,j)*G.y(j)/alpha))/(1.0-exp(v(i,j)/alpha)));
 	}
 }		

void SolveConvectionDiffusion(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme, double alpha)
{

	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();
	
	Vector phi(Nx,Ny);
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);
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



	computeTransientMatrix(A, G, dt);
	InitializeVelCD(u,v,G);
	InitializePhiCD(u,v,phi,G, alpha);
	//Compute Residual and Residual norm
	double R0 = 0.0;
	computeDiffusion(id, phi, G, alpha);
	CDS2(fc_Curr, phi, u, v, G);
	R=id + 1.5 * fc_Curr - 0.5 * fc_Prev;
	applyBC(R,dphi,G);
	R0 = R.L2Norm();
	printf("Initial Residual Norm = %14.12e\n", R0);

	
	int itime = 0; 
	
	while(itime*dt <= tf)
	{




		// record the previous convection vector
		fc_Prev = fc_Curr;
		solveGS(dphi,A,R);

		//Update the solution phi
		phi =  phi + dphi;
		computeDiffusion(id, phi, G, alpha);
		CDS2(fc_Curr, phi, u, v, G);
		R=id+1.5*fc_Curr-0.5*fc_Prev;
		applyBC(R,dphi,G);

		
		//calculate the convection vector at the current time step
		// switch(conScheme){
		// 	case 1:
		// 		FOU(fc_Curr, phi, u, v, G);
		// 		break;
		// 	case 2:
		// 		CDS2(fc_Curr, phi, u, v, G);
		// 		break;
		// 	case 3:
		// 		SOU(fc_Curr, phi, u, v, G);
		// 		break;
		// 	default:
		// 		printf("invalid convection scheme.\n");
		// 		exit(0);
		// }
		// switch(timeScheme){
		// 	case 1:
		// 		eulerExp(phi, fc_Curr, dt);
		// 		break;
		// 	case 2:
		// 		itime < 1 ? eulerExp(phi, fc_Curr, dt) : abs2Exp(phi, fc_Curr, fc_Prev, dt);
		// 		break;
		// 	default:
		// 		printf("invalid time-integration scheme.\n");
		// 		exit(0);
		// }

		// time = (last == true ? tf :(itime+1)*dt);
		// if(itime % 1 == 0 || last)
		// 	printf("Time = %lf, maxphi = %14.12e, minphi = %14.12e\n",time,phi.GetMax(),phi.GetMin());

		// if(last)
		// 	break;

		// ++itime;
		// if(itime*dt > tf){
		// 	dt = tf - (itime-1)*dt;
		// 	last = true;
		// }

		//const double dt2=0.05*2/200;


		//Solve the linear system Adu = -R
		//solveGS(dphi,A,R);






		double R1 = R.L2Norm();

		printf("Time-Step = %d\n",++itime); 
		printf("Residual Norm = %14.12e,\n Residual Norm Ratio (R/R0) = %14.12e\n", R1, R1/R0);
		//Check convergence
		if(dphi.L2Norm() < 1e-8)
		{
			printf("Steady state reached in %d time steps.\n Final time = %lf.\n",itime,itime*dt);
			break;
		}
	}
	char fname[20] = "Parsa.vtk";
	storeVTKStructured(phi, G, fname);
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
		double tf1 = 2*PI;
		double dt1 = 0.05*2/200;
		Grid G1(Nx1,Ny1,xlim1,ylim1);
		SolveConvection(G1, tf1, dt1, conScheme, timeScheme);

	} else if (pbtype == 2){

		std::cout << "2D convection-diffusion equation is solved" << std::endl;
		unsigned long Nx2, Ny2;
		double xlim2[2] = {0, 1};
		double ylim2[2] = {0, 1};
		double alpha = 0.1;
		Nx2 = Ny2 = 17;
		double tf2 = 2*PI;//0.05*2/200;//2*PI;
		double dt2 = 0.05/200;
		Grid G2(Nx2,Ny2,xlim2,ylim2);
		SolveConvectionDiffusion(G2, tf2, dt2, conScheme, timeScheme, alpha);

	} else {
	 	std::cout << "Problem undefined" << std::endl;


	}


}