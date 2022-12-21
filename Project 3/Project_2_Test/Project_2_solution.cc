#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include "Base.h"

#include <iostream>
using namespace std;

double f1(double x, double y);
double f0(double x, double y);
/* add couple of void for matrix and grid */

/*================================================================================================
 * Correct the initialization function for the 2D convection equation f0
 *================================================================================================
 */
double f0(double x, double y)
{
	return 5*exp(-1500*(pow((x-0.25),2) + pow(y,2)));
}
/*================================================================================================*/
void InitializePhi(Vector &phi, const Grid &G)
{
	unsigned long i,j;
	double x, y;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();


    for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++){
			x = G.x(i);
            y = G.y(j); 
			phi(i,j) = f0(x,y);
		}
};
void InitializeVel(Vector &u, Vector &v, const Grid &G)
{
	unsigned long i,j;
	double x,y;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();

	    for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++){
			x = G.x(i);
            y = G.y(j);
			u(i,j)=y;
            v(i,j)=-x;
		}
};
/*================================================================================================
 * Correct the initialization function for the 2D convection diffusion equation f1
 * pay attention to boundary values (use the value from the boundary condition)
 *================================================================================================
 */
double f1(double x, double y)
{
	return exp(-1500*(pow((x-0.5),2) + pow((y-0.5),2)));
};

double fexact(double x, double y, double u, double v, double alpha)
{
	return (5.0*((1-exp(u*x/alpha))/(1-exp(u/alpha)))) + (0.1*((1-exp(y*v/alpha))/(1-exp(v/alpha))));
};
// double f2(double x, double y, double u, double v)
// {
// 	return 5.0*((1-exp((u*x)/alpha))/((1-exp(u/alpha))));  
// };
/*================================================================================================*/
void InitializePhiCD(Vector &phi, const Grid &G)
{
	unsigned long i,j;
	double x,y;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();
    double u = 1;
    double v = 0.05;
	double alpha=0.1; 

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			x = G.x(i);
            y = G.y(j); 
			if(i==0)
			{
				phi(i,j) = 0.1*((1-exp((v*y)/alpha))/((1-exp(v/alpha))));
			}
			else if(i==Nx-1)
			{
				phi(i,j) =5+0.1*((1-exp((v*y)/alpha))/((1-exp(v/alpha))));
			}
            else if(j==0)
            {
                x = G.x(i);
                y = G.y(j+1); 
                phi(i,j) = f1(x,y);
            }
            else if(j==Ny-1)
            {
                x = G.x(i);
                y = G.y(j-1); 
                phi(i,j) = f1(x,y);
            }
			else
			{
				phi(i,j) = f1(x,y);
			}
		}
} ;

void InitializePhiExact(Vector &phi_exact, const Grid &G)
{
	unsigned long i,j;
	double x,y;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();
    double u = 1;
    double v = 0.05;
	double alpha=0.1;

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			x = G.x(i);
            y = G.y(j);
			phi_exact(i,j) = fexact(x,y,u,v,alpha);
			
		}
} ;

void InitializeVelCD(Vector &u, Vector &v, const Grid &G)
{
	unsigned long i,j;
	double x,y;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	    for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++){		
			u(i,j)=1.0;
            v(i,j)=0.05;
		}
};

/*================================================================================================
 * Complete the function for applying convective operator
 * (1) first order upwind (FOU)
 * (2) central difference (CDS2)
 * (3) second order upwind (SOU)
 *================================================================================================
 */
void FOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G) //const vector of u and v?
{
    unsigned long i,j;
	double x,y;
// u,v must be defined
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();

	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++)

        {

        double u_p= 0.5 * (u(i,j) + abs(u(i,j)));
		double u_n= 0.5 * (u(i,j) - abs(u(i,j)));
		double v_p= 0.5 * (v(i,j) + abs(v(i,j)));
		double v_n= 0.5 * (v(i,j) - abs(v(i,j)));

			fc_Curr(i,j) = (u_p/dx)*(phi(i,j)-phi(i-1,j)) + (u_n/dx)*(phi(i+1,j)-phi(i,j)) + (v_p/dy)*(phi(i,j)-phi(i,j-1)) + (v_n/dy)*(phi(i,j+1)-phi(i,j));
		}

};


void CDS2(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
    //is it only for the boundary
    unsigned long i,j;
	double x,y;
// u,v must be defined
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();



	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++)

        {
			fc_Curr(i,j)=(u(i,j)/(2*dx))*(phi(i+1,j)-phi(i-1,j)) + (v(i,j)/(2*dy))*(phi(i,j+1)-phi(i,j-1)); // ALA ERROR
		}
};

void SOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
	unsigned long i,j;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();



for(i = 1; i < Nx-1; i++)
for(j = 1; j < Ny-1; j++)
{

    if ((1<i<Nx-2) && (1<j<Nx-2))
	{


        double u_p= 0.5 * (u(i,j) + abs(u(i,j)));
		double u_m= 0.5 * (u(i,j) - abs(u(i,j)));
		double v_p= 0.5 * (v(i,j)+ abs(v(i,j)));
		double v_m= 0.5 * (v(i,j) - abs(v(i,j)));


		fc_Curr(i,j) = u_p*((3*phi(i,j)-4*phi(i-1,j)+phi(i-2,j))/(2*dx))+u_m*((-3*phi(i,j)+4*phi(i+1,j)-phi(i+2,j))/(2*dx))+v_p*((3*phi(i,j)-4*phi(i,j-1)+phi(i,j-2))/(2*dy))+v_m*((-3*phi(i,j)+4*phi(i,j+1)-phi(i,j+2))/(2*dy));
	}

else

    {fc_Curr(i,j) = -(u(i,j)/(2*dx))*(phi(i+1,j)-phi(i-1,j)) - (v(i,j)/(2*dy))*(phi(i,j+1)-phi(i,j-1)); };


};};

/*================================================================================================
 * Complete the function for time-integration
 * (1) forward Euler (eulerExp)
 * (2) Adams-Bashforth second order (ab2Exp)
 *================================================================================================
 */
void eulerExp(Vector &phi, const Vector &fc_Curr, double &dt)
{
    phi = phi - (fc_Curr*dt);
		
}
void abs2Exp(Vector &phi, const Vector &fc_Curr, const Vector &fc_Prev, double &dt)
{
	phi = phi - ((1.5*fc_Curr*dt)-(0.5*fc_Prev*dt));	
}
void SolveConvection(const Grid &G, const double tf, double dt, const unsigned short scheme);
/*================================================================================================
 * Complete the function for solving convection diffusion equation
 * Variables which might be useful for you has been define inside the function
 *================================================================================================
 */
void SolveConvectionDiffusion(const Grid &G, double tf, double dt, const unsigned short scheme);
/*================================================================================================*/



void SolveConvection(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	Vector phi(Nx,Ny);
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);
	InitializePhi(phi,G);
    char fname1[30] = "Initial_Phi.vtk";
	storeVTKStructured(phi, G, fname1);

	InitializeVel(u,v,G);
    char fname2[30] = "Initial_vel.vtk";
	storeVTKStructured(u, G, fname2);

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

	char fname3[20] = "Phi.vtk";
	storeVTKStructured(phi, G, fname3);


}


/*============================================================================
 * This is the function for constructing the LHS matrix for transient solution
 * Using TRAPEZOIDAL
 * a, b, M(i,j,0), M(i,j,1), M(i,j,2), M(i,j,3), M(i,j,4)
 * You can use dx, dy, dt and operators +,-,*,/ in the calculation
 *============================================================================
*/


void computeTransientMatrix(Matrix &M, const Grid &G, const double &dt)
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

void computeResidual(Vector &R, const Vector &fc_Curr, const Vector &fc_Prev, const Vector &id)
{

    R = - (1.5*fc_Curr) + (0.5*fc_Prev) + id;

}

void computeDiffusion(Vector &id, const Vector &phi, const Grid &G)
{
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();
    double alpha = 0.1;

	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++)
			id(i,j)=alpha*(((phi(i-1,j)-2*phi(i,j)+phi(i+1,j))/(dx*dx))+((phi(i,j-1)-2*phi(i,j)+phi(i,j+1))/(dy*dy)));  // ADD ALPHA HERE

}


void applyBC(Vector &R, Vector &dphi, const Grid &G, const Vector &phi)
{
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();

    double u = 1;
    double v = 0.05;
	double alpha=0.1; 

	const double a = 0.0;

    double dx = G.dx();
	double dy = G.dy();

    const double h1 = 0.1*((-v)/(alpha*(1-exp(v/alpha))));  
    const double h2 = 0.1*((-v*exp(v/alpha))/(alpha*(1-exp(v/alpha))));

	for(i = 0; i < Nx; i++){
		j = 0;

        R(i,j) = h1 - (phi(i,j+1) - phi(i,j))/dy;

		j = Ny-1;
        R(i,j) = h2 - (phi(i,j) - phi(i,j-1))/dy;
	}

	for(j = 0; j < Ny; j++){
		i = 0;
		R(i,j) = dphi(i,j) = a;

		i = Nx-1;
		R(i,j) = dphi(i,j) = a;
	}
}

void SolveConvectionDiffusion(const Grid &G, const double tf, double dt, const unsigned short conScheme, const unsigned short timeScheme)
{

	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	Vector phi(Nx,Ny);
	Vector phi_exact(Nx,Ny);
	Vector u(Nx, Ny);
	Vector v(Nx, Ny);

	InitializePhiCD(phi,G);
	InitializePhiExact(phi_exact,G);
	InitializeVelCD(u,v,G);

    char fname1[30] = "Initial_Phi.vtk";
	storeVTKStructured(phi, G, fname1);

    char fname2[30] = "Initial_vel.vtk";
	storeVTKStructured(v, G, fname2);

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

    //Compute Residual and Residual norm
	double R0 = 0.0;
	computeDiffusion(id,phi,G);
    CDS2(fc_Curr, phi, u, v, G);
    fc_Prev = fc_Curr;
    computeResidual(R,fc_Curr,fc_Prev,id);
	applyBC(R,dphi,G, phi);
	R0 = R.L2Norm();
	printf("Initial Residual Norm = %14.12e\n", R0);


    int itime = 0; double time = 0.0; bool last = false;
	while(time < tf)
	{   
        fc_Prev = fc_Curr;

        solveGS(dphi,A,R);
		phi = phi + dphi;


        // calculate the diffusion vector
        computeDiffusion(id,phi,G);
		// calculate the convection vector at the current time step
		CDS2(fc_Curr, phi, u, v, G);
        // calculate residual
        computeResidual(R,fc_Curr,fc_Prev,id);
        applyBC(R,dphi,G, phi);

        
        double R1 = R.L2Norm();
        
        // if(fmod(itime,500) == 0)
        // {char buffedsr[32]; // Time series
    
    	// snprintf(buffedsr, sizeof(char) * 32, "phi%i.vtk", itime);
		// storeVTKStructured(phi, G, buffedsr);}


		time = (last == true ? tf :(itime+1)*dt);
		if(itime % 1 == 0 || last)
			{printf("Time = %lf, maxphi = %14.12e, minphi = %14.12e\n",time,phi.GetMax(),phi.GetMin());
            printf("Residual Norm = %14.12e,\n Residual Norm Ratio (R/R0) = %14.12e\n", R1, R1/R0);
			printf("Convergence criteria = %14.12e\n", dphi.L2Norm());
            std::cout<<std::endl;}

		if(last)
			break;

		//Check convergence
		if(dphi.L2Norm() < 1e-8){
			printf("Steady state reached in %d time steps.\n Final time = %lf.\n",itime,itime*dt);
			break;
		}

		++itime;
		if(itime*dt > tf){
			dt = tf - (itime-1)*dt;
			last = true;
		}
	}

	char fname3[20] = "Phi.vtk";
	storeVTKStructured(phi, G, fname3);

	char fname4[20] = "Phi_e.vtk";
	storeVTKStructured(phi_exact, G, fname4);

	Vector e(Nx, Ny);
	e = phi - phi_exact;
	printf("L2 error = %14.12e\n", e.L2Norm());


}

int main()
{
	
	unsigned short conScheme = 2;
	unsigned short timeScheme = 1;
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
		Nx2 = Ny2 = 17;
		double tf2 = 50*PI;
		double dt2 = 0.999*1/(Nx2-1);
		Grid G2(Nx2,Ny2,xlim2,ylim2);
		SolveConvectionDiffusion(G2, tf2, dt2, conScheme, timeScheme);

	} else {
	 	std::cout << "Problem undefined" << std::endl;
	}
}