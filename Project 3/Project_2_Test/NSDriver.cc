#include "NSDriver.h"

double test = 0.0;
/*========================================================================================
 * The entire code is written as a class, which takes the input from the NS.inp
 * If you have completed the code correctly, the current set of parameters should give you
 * correct solution. 
 * You can adjust them for testing
 *========================================================================================
 */
/*========================================================================================
 * You need to complete the functions according their name. 
 * You can find the comments above the functions needs to be completed
 * You can find the purposes of the functions above the functions needs to be completed
 *========================================================================================
 */ 
 /*========================================================================================
 * we've updated the linear solver. Now you can solve Ax=b by S.solve(x, A, b);
 *========================================================================================
 */ 


NSParameters::NSParameters():	dt(Config::Instance()->dt()),
								Re(Config::Instance()->Re()),
								nTimeSteps(Config::Instance()->nTimeSteps()),
								alphaP(Config::Instance()->urfPres()),
								alphaU(Config::Instance()->urfMom()),
								convergence(Config::Instance()->conPresMom()),
								maxOuterIterations(Config::Instance()->maxOuterIter())
{}

NSDriver::NSDriver():	G(Config::Instance()->Nx(), Config::Instance()->Ny()),
						Params(NSParameters()),
						S(Solver()), 
						Au(G.Nx(),G.Ny()), Av(G.Nx(),G.Ny()), Ap(G.Nx(),G.Ny()),
						u(G.Nx(),G.Ny()), v(G.Nx(),G.Ny()), p(G.Nx(),G.Ny()),
						u_(u), v_(v), p_(p), Fcu(u), Fcv(v), Fdu(u), Fdv(v),
						gradPx(p), gradPy(p), pc(p), uf(u), vf(v),
						bu(u), bv(u), bp(u),
						bu0(bu), bv0(bv)

{}

NSDriver::~NSDriver()
{}

void NSDriver::initializeVel()
{
	for(size_t i = 0; i < G.Nx(); i++)
		for(size_t j = 0; j < G.Ny(); j++){
			u(i,j) = 0;
			v(i,j) = 0;
		}

	for(size_t i = 0; i < G.Nx(); i++){
		u(i,0) = 0;
		u(i,G.Ny()-1) = 1;
	}
}

void NSDriver::initializePres()
{
	for(size_t i = 0; i < G.Nx(); i++)
		for(size_t j = 0; j < G.Ny(); j++)
			p(i,j) = 0;
}


void NSDriver::Run()
{
	unsigned long itime = 0;
	unsigned long k = 0;
	char fileName[50] = "solution_0.vtk";
	std::cout << std::endl << std::endl << std::endl << Params.Re << std::endl << std::endl << std::endl;
	initializeVel();
	initializePres();

	storeVTKSolution(u,v,p,G,fileName);
	/*========================================================================================
	 * assemble the LHS for the X momentum equation Au (already defined in the sturcture) 
	 *========================================================================================
 	 */ 
	assembleA('U');
	/*========================================================================================
	 * assemble the LHS for the Y momentum equation Av (already defined in the sturcture)
	 *========================================================================================
 	 */
	assembleA('V');

	while(itime < Params.nTimeSteps){
		printf("========================================================\n");
		printf("Time step = %lu, Current Solution time = %f\n",itime+1, (itime+1)*Params.dt);
		printf("========================================================\n");
		// Fcu_ flux of convection term for u at n-1
		Fcu_ = Fcu;
	    /*========================================================================================
	     * given u and v, calculate udu/dx+vdu/dy in conservative form
	     *========================================================================================
 	     */
		computeConvectiveXFlux(Fcu, u, v, G);
	    /*========================================================================================
	     * given u and v, calculate d^2u/dx^2 with central difference scheme
	     *========================================================================================
 	     */
		computeDiffusiveFlux(Fdu, u, G);
	    /*========================================================================================
	     * assemble the residual without the pressure term for X momentum
	     *========================================================================================
 	     */
		// assembleb0('U')
		assembleb0('U', Fcu, Fcu_, Fdu, u, itime);
		// Fcv_ flux of convection term for v at n-1
		Fcv_ = Fcv;
	    /*========================================================================================
	     * given u and v, calculate udv/dx+vdv/dy in conservative form
	     *========================================================================================
 	     */
		computeConvectiveYFlux(Fcv, u, v, G);
	    /*========================================================================================
	     * given u and v, calculate d^2u/dx^2 with central difference scheme
	     *========================================================================================
 	     */
		computeDiffusiveFlux(Fdv, v, G);
	    /*========================================================================================
	     * assemble the residual without the pressure term for Y momentum
	     *========================================================================================
 	     */
		// assembleb0('V');
		assembleb0('V', Fcv, Fcv_, Fdv, v, itime);

		k = 0;

		//Set intermediate flow variables prior to outer iterations
		// Initialize predictors
		u_ = u; v_ = v; p_ = p;

		//Run outer iterations
		while(k < Params.maxOuterIterations){
			printf("----------------------Iteration = %lu-------------------------\n",k);

		/*========================================================================================
	     * given p calculate pressure gradient
	     *========================================================================================
 	     */
			computeGradient(gradPx, gradPy, p_, G);
		/*========================================================================================
	     * predict intermediate u_ and v_ according to guessed p_
	     *========================================================================================
 	     */			
			solveXMomPredictor();
			solveYMomPredictor();
		/*========================================================================================
	     * solve the Poisson's equation for the pressure correction pc
	     * For simplicity, solve it explicitly. No Newton-Raphson iteration is needed
	     * - (1/a)\Delta p = - (du/dx+du/dy+d^4p/dx^4+d^4p/dy^4)
	     * Note that we have handeled the boundary condition for you.
	     * You just need to take care of internal nodes
	     *========================================================================================
 	     */
			solvePresCorrector();
		/*========================================================================================
	     * update the pressure p_ according to pc
	     * update the velocity u_,v_ 
	     * (note that u_ and v_ will be overwritten in the next outer iteration.
	     * u_ and v_ will be passed to the next step when the outer iteration converges.)
	     *========================================================================================
 	     */
			correctFlowVars();
			k++;
			printf("--------------------------------------------------------------\n");
			if(checkConvergence())
				break;
		}

		if(k < Params.maxOuterIterations)
			printf("pressure-velocity coupling converged in %lu iterations.\n",k);
		else
			printf("pressure-velocity coupling limited to %lu iterations.\n",k);

		if((itime+1) % 50 == 0){
			sprintf(fileName,"solution_%lu.vtk",itime+1);
			storeVTKSolution(u,v,p,G,fileName);
		}

		//Update velocity and pressure for next time step
		u = u_; v = v_; p = p_;
		itime++;
	}
}

void NSDriver::solveXMomPredictor()
{
	printf("--------Solving x-Momentum predictor---------\n");


	momResx = bu - Au*u_;
	printf("(Max,Min) predictor velocity = %14.12e, %14.12e\n",u_.GetMax(),u_.GetMin());
}

void NSDriver::solveYMomPredictor()
{
	printf("--------Solving y-Momentum predictor---------\n");


	momResy = bv - Av*v_;
	printf("(Max,Min) predictor velocity = %14.12e, %14.12e\n",v_.GetMax(),v_.GetMin());
}

void NSDriver::solvePresCorrector()
{
	printf("--------Solving Pressure Corrector---------\n");

	assembleAp();
	assemblebp();
	S.solve(pc, Ap, bp);
	printf("(Max,Min) pressure correction = %14.12e, %14.12e\n",pc.GetMax(),pc.GetMin());
}

void NSDriver::correctFlowVars()
{
	double dx = G.dx(), dy = G.dy();
	const double a = (1.0/Params.dt + 1.0/(Params.Re*dx*dx) + 1.0/(Params.Re*dy*dy));

	p_ = p_ + Params.alphaP*pc;
	for(size_t i = 1; i < G.Nx()-1; i++)
		for(size_t j = 1; j < G.Ny()-1; j++){
			// complete the following two lines
			u_(i,j) = u_(i,j) - 0;
			v_(i,j) = v_(i,j) - 0;
		}
}

bool NSDriver::checkConvergence()
{
	double presCorrNorm = pc.L2Norm();
	double momResNorm[2] = {momResx.L2Norm(), momResy.L2Norm()};
	printf("Pressure correction norm = %14.12e\n",presCorrNorm);
	printf("Momentum residuals = %14.12e, %14.12e\n",momResNorm[0],momResNorm[1]);

	if(MAX(presCorrNorm,MAX(momResNorm[0],momResNorm[1])) < Params.convergence)
		return true;
	return false;
}

void NSDriver::assembleb0(char flowVar,const Vector &Fc, const Vector &Fc_,const Vector &Fd, const Vector &vel, unsigned long itime)
{
	switch(flowVar)
	{
		case 'U':
		// Calculate non-convective terms first
		bu0 = (1 / Params.dt) * vel + (0.5 / Params.Re) * Fd;
		// Foward Euler if first timestep, Adams Bashforth otherwise 
		if (itime < 1)
		{
			bu0 = bu0 - Fc;
		}
		else
		{
			bu0 = bu0 - (1.5 * Fc - 0.5 * Fc_);
		}
		break;
		case 'V':
		// Calculate non-convective terms first
		bu0 = (1 / Params.dt) * vel + (0.5 / Params.Re) * Fd;
		// Foward Euler if first timestep, Adams Bashforth otherwise 
		if (itime < 1)
		{
			bu0 = bu0 - Fc;
		}
		else
		{
			bu0 = bu0 - (1.5 * Fc - 0.5 * Fc_);
		}
		break;
	}
}

void NSDriver::assembleA(char flowVar)
{
	double Vc = G.dx()*G.dy();
	double dx, dy;
	dx = G.dx(), dy = G.dy();

	switch(flowVar)
	{
		case 'U':
		for(size_t i = 1; i < G.Nx()-1; i++)
			for(size_t j = 1; j < G.Ny()-1; j++)
			{
				Au(i,j,0) = 0;
				// Added code below here to try
				Au(i, j, 1) = 0;
				Au(i, j, 2) = ((1 / Params.dt) - (0.5 / Params.Re));
				Au(i, j, 3) = 0;
				Au(i, j, 4) = 0;
			}	
		break;

		case 'V':
		for(size_t i = 1; i < G.Nx()-1; i++)
			for(size_t j = 1; j < G.Ny()-1; j++)
			{
				Av(i,j,0) = 0;
				// added code below here to try
				Av(i, j, 1) = 0;
				Av(i, j, 2) = ((1 / Params.dt) - (0.5 / Params.Re));
				Av(i, j, 3) = 0;
				Av(i, j, 4) = 0;
			}	
		break;
	}
}


void NSDriver::assembleAp()
{
	//TODO ADD COMPUTE TRANSIENT MATRIX AND COMPUTE RESIDUAL FROM PROJECT 2 TO HERE AND ASSEMBLEBP
	double dx, dy;
	dx = G.dx(), dy = G.dy();
	size_t Nx, Ny;
	Nx = G.Nx(); Ny = G.Ny();
	const double a = (1.0/Params.dt + 1.0/(Params.Re*dx*dx) + 1.0/(Params.Re*dy*dy));
	double fct = 0;
	for (size_t i = 2; i < Nx-2; i++){
		for (size_t j = 2; j < Ny-2; j++){
		// Ap for internal nodes
		// Rewrite the entire function for testing Project 1
		// Copy the right value to the following five lines.
		// For solving NS equation, changing the rest of the code is not recommended.
		
			Ap(1,j,0) = test;
			Ap(1,j,1) = test;
			Ap(1,j,2) = test;
			Ap(1,j,3) = test;
			Ap(1,j,4) = test;
		}
	}
	for (size_t j = 2; j < Ny-2; j++){
			Ap(1,j,0) = - fct*1.0/(a*dx*dx);
			Ap(1,j,1) = - 1.0/(a*dy*dy);
			Ap(1,j,2) =   (1+fct)/(a*dx*dx)+2.0/(a*dy*dy);
			Ap(1,j,3) = - 1.0/(a*dy*dy);
			Ap(1,j,4) = - 1.0/(a*dx*dx);

			Ap(Nx-2,j,0) = - 1.0/(a*dx*dx);
			Ap(Nx-2,j,1) = - 1.0/(a*dy*dy);
			Ap(Nx-2,j,2) =   (1+fct)/(a*dx*dx)+2.0/(a*dy*dy);
			Ap(Nx-2,j,3) = - 1.0/(a*dy*dy);
			Ap(Nx-2,j,4) = - fct*1.0/(a*dx*dx);
	}

		for (size_t i = 2; i < Nx-2; i++){
			Ap(i,1,0) = - 1.0/(a*dx*dx);
			Ap(i,1,1) = - fct*1.0/(a*dy*dy);
			Ap(i,1,2) =   2.0/(a*dx*dx)+(1+fct)/(a*dy*dy);
			Ap(i,1,3) = - 1.0/(a*dy*dy);
			Ap(i,1,4) = - 1.0/(a*dx*dx);

			Ap(i,Ny-2,0) = - 1.0/(a*dx*dx);
			Ap(i,Ny-2,1) = - 1.0/(a*dy*dy);
			Ap(i,Ny-2,2) =   2.0/(a*dx*dx)+(1+fct)/(a*dy*dy);
			Ap(i,Ny-2,3) = - fct*1.0/(a*dy*dy) ;
			Ap(i,Ny-2,4) = - 1.0/(a*dx*dx);
	}

			Ap(1,1,0) = - fct*1.0/(a*dx*dx) ;
			Ap(1,1,1) = - fct*1.0/(a*dy*dy) ;
			Ap(1,1,2) =   (1+fct)/(a*dx*dx)+(1+fct)/(a*dy*dy);
			Ap(1,1,3) = - 1.0/(a*dy*dy);
			Ap(1,1,4) = - 1.0/(a*dx*dx);

			Ap(1,Ny-2,0) = - fct*1.0/(a*dx*dx);
			Ap(1,Ny-2,1) = - 1.0/(a*dy*dy);
			Ap(1,Ny-2,2) =   (1+fct)/(a*dx*dx)+(1+fct)/(a*dy*dy);
			Ap(1,Ny-2,3) = - fct*1.0/(a*dy*dy);
			Ap(1,Ny-2,4) = - 1.0/(a*dx*dx);

			Ap(Nx-2,1,0) = - 1.0/(a*dx*dx);
			Ap(Nx-2,1,1) = - fct*1.0/(a*dy*dy);
			Ap(Nx-2,1,2) =   (1+fct)/(a*dx*dx)+(1+fct)/(a*dy*dy);
			Ap(Nx-2,1,3) = - 1.0/(a*dy*dy);
			Ap(Nx-2,1,4) = - fct*1.0/(a*dx*dx);

			Ap(Nx-2,Ny-2,0) = - 1.0/(a*dx*dx);
			Ap(Nx-2,Ny-2,1) = - 1.0/(a*dy*dy);
			Ap(Nx-2,Ny-2,2) =   (1+fct)/(a*dx*dx)+(1+fct)/(a*dy*dy);
			Ap(Nx-2,Ny-2,3) = - fct*1.0/(a*dy*dy);
			Ap(Nx-2,Ny-2,4) = - fct*1.0/(a*dx*dx);

			for (size_t j = 1; j < Ny-1; j++){
				Ap(0,j,0) =   0;
				Ap(0,j,1) =   0;
				Ap(0,j,2) =   1.0/(a*dx*dy);
				Ap(0,j,3) =   0;
				Ap(0,j,4) = - 1.0/(a*dx*dy);

				Ap(Nx-1,j,0) = - 1.0/(a*dx*dy);
				Ap(Nx-1,j,1) =   0;
				Ap(Nx-1,j,2) =   1.0/(a*dx*dy);
				Ap(Nx-1,j,3) =   0;
				Ap(Nx-1,j,4) =   0;
			}

			for (size_t i = 1; i < Nx-1; i++){
				Ap(i,0,0) =   0;
				Ap(i,0,1) =   0;
				Ap(i,0,2) =   1.0/(a*dx*dy);
				Ap(i,0,3) = - 1.0/(a*dx*dy);
				Ap(i,0,4) =   0;

				Ap(i,Ny-1,0) =   0;
				Ap(i,Ny-1,1) = - 1.0/(a*dx*dy);
				Ap(i,Ny-1,2) =   1.0/(a*dx*dy);
				Ap(i,Ny-1,3) =   0;
				Ap(i,Ny-1,4) =   0;
			}

			Ap(0,0,0) = 0;
			Ap(0,0,1) = 0;
			Ap(0,0,2) = 1.0/(a*dx*dy);
			Ap(0,0,3) = 0;//-1.0/(a*dx*dy);
			Ap(0,0,4) = 0;

			Ap(Nx-1,0,0) = 0;
			Ap(Nx-1,0,1) = 0;
			Ap(Nx-1,0,2) = 1.0/(a*dx*dy);
			Ap(Nx-1,0,3) = 0;
			Ap(Nx-1,0,4) = 0;

			Ap(0,Ny-1,0) = 0;
			Ap(0,Ny-1,1) = 0;//-1.0/(a*dx*dy);
			Ap(0,Ny-1,2) = 1.0/(a*dx*dy);
			Ap(0,Ny-1,3) = 0;
			Ap(0,Ny-1,4) = 0;
			
			Ap(Nx-1,Ny-1,0) = 0;
			Ap(Nx-1,Ny-1,1) = 0;//-1.0/(a*dx*dy);
			Ap(Nx-1,Ny-1,2) = 1.0/(a*dx*dy);
			Ap(Nx-1,Ny-1,3) = 0;
			Ap(Nx-1,Ny-1,4) = 0;

			Ap(Nx-2,1,0) = 0;
			Ap(Nx-2,1,1) = 0;//-1.0/(a*dx*dy);
			Ap(Nx-2,1,2) = 1.0/(a*dx*dy);
			Ap(Nx-2,1,3) = 0;
			Ap(Nx-2,1,4) = 0;

}


void NSDriver::assemblebp()
{
	double dx, dy;
	dx = G.dx(), dy = G.dy();
	size_t Nx, Ny;
	Nx = G.Nx(); Ny = G.Ny();
	const double a = (1.0/ Params.dt + 1.0/(Params.Re*dx*dx) + 1.0/(Params.Re*dy*dy));
	// extrapolated pressure at east, west, north and south
	double ee, ww, nn, ss;
	double et, wt, nt, st;
	for (size_t i = 1; i < Nx-1; i++){
		for (size_t j = 1; j < Ny-1; j++){
			// tilde = 0
			// et = p_(Nx-2,j);
			// wt = p_(1,j);
			// st = p_(i,1);
			// nt = p_(i,Ny-2);

			// 0th order
			// et = p_(Nx-1,j);
			// wt = p_(0,j);
			// st = p_(i,0);
			// nt = p_(i,Ny-1);

			// 1st order
			// et = 2.0 * p_(Nx-1,j) - p_(Nx-2,j);
			// wt = 2.0 * p_(0,j)    - p_(1,j);
			// st = 2.0 * p_(i,0)    - p_(i,1);
			// nt = 2.0 * p_(i,Ny-1) - p_(i,Ny-2);

			// 2nd order
			wt = 2.5 * p_(0,j)    -2  * p_(1,j)    + 0.5*p_(2,j);
			et = 2.5 * p_(Nx-1,j) - 2 * p_(Nx-2,j) + 0.5*p_(Nx-3,j);
			nt = 2.5 * p_(i,Ny-1) - 2 * p_(i,Ny-2) + 0.5*p_(i,Ny-3);
			st = 2.5 * p_(i,0)    - 2 * p_(i,1)    + 0.5*p_(i,2);
			
			if (i == 1 && j != 1 && j != Ny-2){
				ww = wt ;
				ee = p_(i+2,j);
				nn = p_(i,j+2);
				ss = p_(i,j-2);
			}
			else if (i == Nx-2 && j != 1 && j != Ny-2){ 
				ww = p_(i-2,j);
				ee = et ;
				nn = p_(i,j+2);
				ss = p_(i,j-2);
			}
			else if (j == 1 && i != 1 && i != Nx-2){
				ww = p_(i-2,j);
				ee = p_(i+2,j);
				nn = p_(i,j+2);
				ss = st ;
			}
			else if (j == Ny-2 && i != 1 && i != Nx-2){
				ww = p_(i-2,j);
				ee = p_(i+2,j);
				nn = nt;
				ss = p_(i,j-2);
			}
			else if (i == 1 && j == 1){
				ww = wt; 
				ee = p_(i+2,j);
				nn = p_(i,j+2);
				ss = st;
			}	
			else if (i == 1 && j == Ny-2){
				ww = wt; 
				ee = p_(i+2,j);
				nn = nt;
				ss = p_(i,j-2);
			}
			else if (i == Nx-2 && j == 1){
				ww = p_(i-2,j);
				ee = et;
				nn = p_(i,j+2);
				ss = st;
			}
			else if (i == Nx-2 && j == Ny-2){
				ww = p_(i-2,j);
				ee = et;
				nn = nt;
				ss = p_(i,j-2);
			} else {
				ww = p_(i-2,j);
				ee = p_(i+2,j);
				nn = p_(i,j+2);
				ss = p_(i,j-2);
			}
			// bp for internal nodes
			// You just need to assign value for bp in the next line.
			// Rewrite the entire function for testing Project 1
			// For solving NS equation, changing the rest of the code is not recommended.
			// use ww for p(i-2,j), use ee for p(i+2,j)
			// use nn for p(i,j+2), use ss for p(i,j-2)  
			bp(i,j) = test;

		}
	}
	for (size_t i = 0; i < Nx; i++){
		bp(i,0) = 0;
		bp(i,Ny-1) = 0;
	}
		
	for (size_t j = 0; j < Ny; j++){
		bp(0,j) = 0;
		bp(Nx-1,j) = 0;	
	}
		
}


int main()
{
	NSDriver NS;
	NS.Run();	

}