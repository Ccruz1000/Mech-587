#include <iostream>
#include "Base.h"

// Input Variables
double conv_u = 1.0; // Convective u velocity cx
double conv_v = 0.05; // Convective v velocit cv
double Re = 10;

// // Function for the initial value from project 2
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
void initvelexact(const Grid &G)
{
	unsigned long i, j;
	double x, y;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	Vector vel_exact(Nx, Ny);
	const double dx = G.dx();
	const double dy = G.dy();
	double alpha = 1 / Re;
	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			x = G.x(i);
            y = G.y(j);
			vel_exact(i,j) = fexact(x, y, conv_u, conv_v, alpha);
			
		}
	char fname[20] = "vel_exact.vtk";
	storeVTKStructured(vel_exact, G, fname);
}

// // Initialize the convective velocity Cx and Cy for part 1
void InitVelConv(Vector &u, Vector &v, const Grid &G)
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
	char uname[20] = "conv_velu_exact.vtk";
	char vname[20] = "conv_velv_exact.vtk";
	storeVTKStructured(u, G, uname);
	storeVTKStructured(v, G, vname);
}

// Here we are going to use the velocity vector as vel from project 2. This is because
// in the momentum equation we are trying to solve for the velocity
void InitVel(Vector &vel, const Grid &G)
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
				vel(i,j) = 0.1*((1-exp((conv_v*y)/alpha))/((1-exp(conv_v/alpha))));
			}
			else if(i==Nx-1)
			{
				vel(i,j) =5+0.1*((1-exp((conv_v*y)/alpha))/((1-exp(conv_v/alpha))));
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

// // Calcconv_ulate first derivative using central difference method 
// void CDS2(Vector &Flux_Conv_Curr, const Vector &vel, const Vector &uconv, const Vector &vconv, const Grid &G)
// {

// }

void SolveConvectionDiffusion(const Grid G)
{
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();

	Vector uconv(Nx, Ny);
	Vector vconv(Nx, Ny);
	initvelexact(G);
	InitVelConv(uconv, vconv, G);
}

int main()
{
	unsigned long Nx, Ny;
	Nx = Ny = 17;
	double xlim[2] = {0, 1};
	double ylim[2] = {0, 1};
	Grid G(Nx, Ny, xlim, ylim);
	SolveConvectionDiffusion(G);

	return 0;
}