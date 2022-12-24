// Include Statements
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include <iostream>
#include "Base.h"

// Input variables
double dt = 0.1;
unsigned long Re = 2;
unsigned long Nstep = 33;
unsigned long Nx = Nstep;
unsigned long Ny = Nstep;
unsigned long nTimeSteps = 500;

// Initialize for project 1
double f(double x, double y)
{
	return pow(x,4) + pow(y,4)- 6*pow(x,2)*pow(y,2);
}
double f2_proj1(double x, double y)
{
	return exp(-100*(pow(x,2) + pow(y,2)));
}


// Initial velocity for project 2
double f1(double x, double y)
{
	return exp(-1500*(pow((x-0.5),2) + pow((y-0.5),2)));
}

// Initialize the pressure for p=x and p=y
void InitializePressure(Vector &px, Vector &py, const Grid &G)
{
	unsigned long i, j;
	double x, y;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();

	for(i = 0; i < Nx; i++)
		for(j = 0; j < Ny; j++)
		{
			x = G.x(i);
			y = G.y(j);
			px(i, j) = x; // For P = x test
			py(i, j) = y; // For P = y test
		}
}

// Calculate the gradient of the pressure for part 1
void CalcGradPressure(Vector &gradpx, Vector &gradpy, const Vector &px,const Vector &py, const Grid &G)
{
	unsigned long i, j;
	const unsigned long Nx = G.Nx();
	const unsigned long Ny = G.Ny();
	const double dx = G.dx();
	const double dy = G.dy();
	// Interior points only
	for(i = 1; i < Nx - 1; i++)
		for(j = 1; j < Ny - 1; j++)
		{
			gradpx(i, j) = (px(i + 1, j) - px(i - 1, j)) / (2 * dx);
			gradpy(i, j) = (py(i, j + 1) - py(i, j - 1)) / (2 * dy); 
		}
}

// function for the exact solution for project 2
double fexact(double x, double y, double u, double v, double alpha)
{
	return (5.0*((1-exp(u*x/alpha))/(1-exp(u/alpha)))) + (0.1*((1-exp(y*v/alpha))/(1-exp(v/alpha))));
};

// Initialize for project 2
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

// Initialize the exact solution for project 2
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

// Initialize convective velocity for project 2
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


void FOU(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G) //const vector of u and v?
{
    unsigned long i,j;
	double x,y;
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

// Central difference for project 2
void CDS2(Vector &fc_Curr, const Vector &phi, const Vector &u, const Vector &v, const Grid &G)
{
    unsigned long i,j;
	double x,y;
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
    	{
    		fc_Curr(i,j) = -(u(i,j)/(2*dx))*(phi(i+1,j)-phi(i-1,j)) - (v(i,j)/(2*dy))*(phi(i,j+1)-phi(i,j-1)); };
		};
};

// Compute matrix to solve for project 1
void computeMatrix(Matrix &M, const Grid &G)
{
	unsigned long i,j;
	const double dx = G.dx();
	const double dy = G.dy();

	for(i = 1; i < G.Nx()-1; i++)
		for(j = 1; j < G.Ny()-1; j++){
			M(i,j,0) = 1.0/(dx*dx);
			M(i,j,1) = 1.0/(dy*dy);
			M(i,j,3) = 1.0/(dy*dy);
			M(i,j,4) = 1.0/(dx*dx);
			M(i,j,2) = -(M(i,j,0) + M(i,j,1) + M(i,j,3) + M(i,j,4));
		}

	for(i = 0; i < G.Nx(); i++)
		for(int t = 0; t < 5; t++){
			M(i,0,t) = (t == 2 ? 1.0 : 0.0);
			M(i,G.Ny()-1, t) = (t == 2 ? 1.0 : 0.0);
		}

	for(j = 0; j < G.Ny(); j++)
		for(int t = 0; t < 5; t++){
			M(0,j,t) = (t == 2 ? 1.0 : 0.0);
			M(G.Nx()-1,j,t) = (t == 2 ? 1.0 : 0.0);
		}
}

// Compute matrix for project 2
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

// Compute residual for project 2
void computeResidual(Vector &R, const Vector &fc_Curr, const Vector &fc_Prev, const Vector &id)
{

    R = - (1.5*fc_Curr) + (0.5*fc_Prev) + id;

}

// computeDiffusion for project 1
void computeDiffusion1(Vector &R, const Vector &u, const Grid &G)
{
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();

	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++)
			R(i,j) = ((u(i+1,j) - 2*u(i,j) + u(i-1,j))/(dx*dx) + (u(i,j+1) - 2*u(i,j) + u(i,j-1))/(dy*dy));

}

// Compute diffusion for project 2
void computeDiffusion2(Vector &id, const Vector &phi, const Grid &G)
{
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();
	double dx = G.dx();
	double dy = G.dy();
    double alpha=0.1;

	for(i = 1; i < Nx-1; i++)
		for(j = 1; j < Ny-1; j++)
			id(i,j)=alpha*(((phi(i-1,j)-2*phi(i,j)+phi(i+1,j))/(dx*dx))+((phi(i,j-1)-2*phi(i,j)+phi(i,j+1))/(dy*dy)));  // ADD ALPHA HERE

}

// Initialize vector for project 1 solution
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
			u(i,j) = f2_proj1(x,y);
		}
}

// Apply boundary conditions for project 1
void applyBC1(Vector &R, Vector &du, const Grid &G)
{
	unsigned long i,j;
	unsigned long Nx = G.Nx();
	unsigned long Ny = G.Ny();

	for(i = 0; i < Nx; i++){
		j = 0;
		R(i,j) = du(i,j) = 0.0;

		j = Ny-1;
		R(i,j) = du(i,j) = 0.0;
	}

	for(j = 0; j < Ny; j++){
		i = 0;
		R(i,j) = du(i,j) = 0.0;

		i = Nx-1;
		R(i,j) = du(i,j) = 0.0;
	}
}


// Apply Boundary conditions for project 2
void applyBC2(Vector &R, Vector &dphi, const Grid &G, const Vector &phi)
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

// Solve steady laplace equation for project 1
const Vector solveSteadyLaplace(const Grid &G, Vector &u)
{
	unsigned long Nx, Ny;
	Nx = G.Nx(), Ny = G.Ny();

	Matrix A(Nx,Ny);
	Vector R(Nx,Ny);
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
	computeDiffusion1(R, u, G);
	applyBC1(R,du,G);
	R = -R;
	R0 = R.L2Norm();
	// char vecName[20] = "SteadyR0.dat";
	// R.storeV(vecName);
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
		computeDiffusion1(R, u, G);
		applyBC1(R, du, G);
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

void SolveConvectionDiffusion(const Grid &G, const double tf, double dt)
{
	const size_t Nx = G.Nx();
	const size_t Ny = G.Ny();

	// Solve and plot exact solution
	Vector phi_exact(Nx,Ny);
	InitializePhiExact(phi_exact, G);
	char fname4[20] = "Phi_exact.vtk";
	storeVTKStructured(phi_exact, G, fname4);

	// Solve and plot initial variables for initial conditions
	Vector u_conv(Nx, Ny);
	Vector v_conv(Nx, Ny);
	Vector px(Nx, Ny); // Pressure = x
	Vector gradpx(Nx, Ny); // Gradient of pressure w.r.t X
	Vector py(Nx, Ny); // Pressure = Y
	Vector gradpy(Nx, Ny); // Gradient of pressure w.r.t. y
	Vector p_sol(Nx, Ny); // Vector to test using project 1
	solveSteadyLaplace(G, p_sol);
	InitializeVelCD(u_conv,v_conv,G); // Initialize Convective velocity
	InitializePressure(px, py, G); // Initialize pressure for part 1
	CalcGradPressure(gradpx, gradpy, px, py, G); // Calculate Pressure Gradient
	// Save Pressure Gradients for problem 1
	char pxvtkfile[30] = "P=X.vtk";
	storeVTKStructured(px, G, pxvtkfile);
	char pyvtkfile[30] = "P=Y.vtk";
	storeVTKStructured(py, G, pyvtkfile);
	char gradpxvtkfile[30] = "gradxP=X.vtk";
	storeVTKStructured(gradpx, G, gradpxvtkfile);
	char gradpyvtkfile[30] = "GradyP=Y.vtk";
	storeVTKStructured(gradpy, G, gradpyvtkfile);
	char fname1[30] = "Init_Conv.vtk";
	storeVTKSolution(u_conv, v_conv, p_sol, G, fname1);

	/**********************************************
	 Declare variables and begin solving for 
	 X-momentum equation 
	 **********************************************/
	Vector VelX(Nx,Ny); // X-velocity 
	InitializePhiCD(VelX,G); // Initialize X velocity
	// variables needed for Newton-Raphson iteration in X
	Matrix Ax(Nx,Ny); // LHS matrix Ax
	Vector Rx(Nx,Ny); // Residual Rx
	Vector dVelX(Nx,Ny); // Increment dphi
	// variables needed for calculating the residual
	Vector fc_Curr_X(Nx,Ny); // Forward convection
	Vector fc_Prev_X(Nx,Ny); // Previous Convection
	Vector diffuseX(Nx,Ny); // Diffusive Residual 
	// Solve First Iteration
    computeTransientMatrix(Ax, G, dt);
	double R0X = 0.0; // X Residual Norm
	computeDiffusion2(diffuseX,VelX,G);
    CDS2(fc_Curr_X, VelX, u_conv, v_conv, G);
    fc_Prev_X = fc_Curr_X;
    computeResidual(Rx,fc_Curr_X,fc_Prev_X,diffuseX);
	applyBC2(Rx,dVelX,G, VelX);
	R0X = Rx.L2Norm();
	printf("Initial X Residual Norm = %14.12e\n", R0X);

	/**********************************************
	 Declare variables and begin solving for 
	 Y-momentum equation 
	 **********************************************/
	Vector VelY(Nx,Ny); // Y-velocity 
	InitializePhiCD(VelY,G); // Initialize Y velocity
	// variables needed for Newton-Raphson iteration in Y
	Matrix Ay(Nx,Ny); // LHS matrix Ay
	Vector Ry(Nx,Ny); // Residual Ry
	Vector dVelY(Nx,Ny); // Increment dphi
	// variables needed for calculating the residual
	Vector fc_Curr_Y(Nx,Ny); // Forward convection
	Vector fc_Prev_Y(Nx,Ny); // Previous Convection
	Vector diffuseY(Nx,Ny); // Diffusive Residual 
	// Solve First Iteration
    computeTransientMatrix(Ay, G, dt);
	double R0Y = 0.0; // X Residual Norm
	computeDiffusion2(diffuseY,VelY,G);
    CDS2(fc_Curr_Y, VelY, u_conv, v_conv, G);
    fc_Prev_Y = fc_Curr_Y;
    computeResidual(Ry,fc_Curr_Y,fc_Prev_Y,diffuseY);
	applyBC2(Ry,dVelY,G, VelY);
	R0Y = Ry.L2Norm();
	printf("Initial Y Residual Norm = %14.12e\n", R0Y);

	char fileName[50] = "solution_0.vtk";
	storeVTKSolution(VelX, VelY, p_sol, G, fileName);


    unsigned long itime = 0; double time = 0.0; bool last = false;
	while(time < tf)
	{   
		// Update convection terms
        fc_Prev_X = fc_Curr_X;
        fc_Prev_Y = fc_Curr_Y;

        // Solve Au=b for X and Y Momentum
        solveGS(dVelX,Ax,Rx);
        solveGS(dVelY,Ay,Ry);

        // Increment solution
		VelX = VelX + dVelX;
		VelY = VelY + dVelY;

        // calculate the diffusion vector for X and Y momentum
        computeDiffusion2(diffuseX,VelX,G);
        computeDiffusion2(diffuseY,VelY,G);

		// calculate the convection vector at the current time step
		CDS2(fc_Curr_X, VelX, u_conv, v_conv, G);
		CDS2(fc_Curr_Y, VelY, u_conv, v_conv, G);

        // Calculate Residual
        computeResidual(Rx,fc_Curr_X,fc_Prev_X,diffuseX);
        computeResidual(Ry,fc_Curr_Y,fc_Prev_Y,diffuseY);

        // Apply Boundary Conditions 
        applyBC2(Rx,dVelX,G, VelX);
        applyBC2(Ry,dVelY, G, VelY);

        // Calculate Residual Norm 
        double R1X = Rx.L2Norm();
        double R1Y = Ry.L2Norm();
        
        // Check if final time step is reached
		time = (last == true ? tf :(itime+1)*dt);
		if(itime % 1 == 0 || last)
		{
			std::cout<<std::endl;
			printf("==================================================================================================\n");
			printf("Timestep = %ld, Time = %lf, X-Residual = %14.12e, Y-Residual = %14.12e\n", itime, time, 
					dVelX.L2Norm(), dVelY.L2Norm());
			printf("Max X-Velocity = %14.12e, Min X-Velocity = %14.12e,\n", 
					VelX.GetMax(),VelX.GetMin());
			printf("X Residual Norm = %14.12e,\nX Residual Norm Ratio (R/R0X) = %14.12e\n", R1X, R1Y/R0X);
			printf("Max Y-Velocity = %14.12e, Min Y-Velocity = %14.12e\n", VelY.GetMax(),VelY.GetMin());
			printf("Y Residual Norm = %14.12e,\nY Residual Norm Ratio (R/R0X) = %14.12e\n", R1Y, R1Y/R0Y);
			printf("==================================================================================================\n");
            std::cout<<std::endl;
        }

		if(last)
			break;

		if((itime + 1) % 10 == 0)
		{
			sprintf(fileName, "solution_%lu.vtk", itime + 1);
			storeVTKSolution(VelX, VelY, p_sol, G, fileName);
		}

		//Check convergence
		if(dVelX.L2Norm() < 1e-8 && dVelY.L2Norm()){
			printf("Steady state reached in %lu time steps.\n Final time = %lf.\n",itime,itime*dt);
			break;
		}

		++itime;
		if(itime*dt > tf){
			dt = tf - (itime-1)*dt;
			last = true;
		}
	}

	Vector e(Nx, Ny);
	e = VelX - phi_exact;
	printf("L2 error = %14.12e\n", e.L2Norm());
}

int main()
{
	double xlim[2] = {0, 1};
	double ylim[2] = {0, 1};
	Grid G2(Nx, Ny, xlim, ylim);
	double tf = nTimeSteps * dt;
	SolveConvectionDiffusion(G2, tf, dt);
}
