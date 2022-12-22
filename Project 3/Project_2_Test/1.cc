#include "NSDriver.h"

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
	for (size_t i = 0; i < G.Nx(); i++){
		for (size_t j = 0; j < G.Ny(); j++){
			p(i, j) = 0;
		}
	}
	p(G.Nx() - 1, G.Ny() - 1) = 1;
}


void NSDriver::Run()
{
	// Input Parameters
	double grid_size = 10;   // For uniform grid dx = dy
	double Nx = grid_size;
	double Ny = grid_size;

	double domain_length = 1;
	double Re = 100;
	double nu = 1 / Re;
	double dt = Params.dt;

	double dx = domain_length / (Nx - 1);
	double dy = domain_length / (Ny - 1);
	double h = domain_length / (Nx - 1);   // For uniform grid dx = dy

	unsigned long i, j;

	// Under relaxation factors:
	double alpha = 0.8;   // for velocity
	double alpha_p = 0.8;   // for pressure
	
	// Initialize variables // test

	Vector u_final(Nx - 1, Ny - 1);
	Vector v_final(Nx - 1, Ny - 1);
	Vector p_final(Nx - 1, Ny - 1);
	
	u_final(Nx - 1, Ny - 1) = 0;
	v_final(Nx - 1, Ny - 1) = 0;
	p_final(Nx - 1, Ny - 1) = 1;

	// Top BC
	for (j = 0; j < Ny; j++) {
		u_final(0, j) = 1;
	}
	
	// Staggered variables
	
	Vector u(Nx, Ny);
	Vector u_star(Nx, Ny);
	Vector d_e(Nx, Ny);
	Vector v(Nx, Ny);
	Vector v_star(Nx, Ny);
	Vector d_n(Nx, Ny);
	Vector p(Nx, Ny);
	Vector p_star(Nx, Ny);
	Vector pc(Nx, Ny);
	Vector b(Nx, Ny);
	
	u(Nx, Ny) = 0;
	u_star(Nx, Ny) = 0;
	d_e(Nx, Ny) = 0;
	v(Nx, Ny) = 0;
	v_star(Nx, Ny) = 0;
	d_n(Nx, Ny) = 0;
	p(Nx, Ny) = 1;
	p_star(Nx, Ny) = 1;
	pc(Nx, Ny) = 0;
	b(Nx, Ny) = 0;

	// Top BC
	for (j = 0; j < Ny; j++) {
		u(0, j) = 2;   // Double because staggered grid (2,0 method)
	}

	Vector u_new(Nx, Ny);
	Vector v_new(Nx, Ny);
	Vector p_new(Nx, Ny);

	u_new(Nx, Ny) = 0;
	v_new(Nx, Ny) = 0;
	p_new(Nx, Ny) = 1;

	// Top BC
	for (j = 0; j < Ny; j++) {
		u_new(0, j) = 2;   // Double because staggered grid (2,0 method)
	}

	// Solve the governing equations
	double error = 1;
	double iterations = 0;
	double error_req = 1e-7;   // final required residual error

	double u_E, u_W, v_N, v_S;
	double a_E, a_W, a_N, a_S;
	double a_e, A_e, a_n, A_n;
	double a_P;


	
	while (error > error_req) {
		
		// x-momentum equation (Interior nodes)
		for (i = 1; i < Nx; i++) {
			for (j = 1; j < Ny - 1; j++) {
				u_E = 0.5 * (u(i, j) + u(i, j + 1));
				u_W = 0.5 * (u(i, j) + u(i, j - 1));
				v_N = 0.5 * (v(i - 1, j) + v(i - 1, j + 1));
				v_S = 0.5 * (v(i, j) + v(i, j + 1));

				a_E = -0.5 * u_E * h + nu;
				a_W = 0.5 * u_W * h + nu;
				a_N = -0.5 * v_N * h + nu;
				a_S = 0.5 * v_S * h + nu;

				a_e = 0.5 * u_E * h - 0.5 * u_W * h + 0.5 * v_N * h - 0.5 * v_S * h + 4 * nu;
				A_e = -h;
				
				d_e(i, j) = A_e / a_e;

				u_star(i, j) = (a_E * u(i, j + 1) + a_W * u(i, j - 1) + a_N * u(i - 1, j) + a_S * u(i + 1, j)) / a_e + d_e(i, j) * (p(i, j + 1) - p(i, j));
			}
		}

		// x-momentum BC's
		for (j = 0; j < Ny; j++) {
			u_star(0, j) = 2 - u_star(1, j);
			u_star(Nx, j) = -u_star(Nx - 1, j);
		}
		
		for (i = 1; i < Nx; i++) {
			u_star(i, 0) = 0;
			u_star(i, Ny - 1) = 0;
		}

		// y-momentum equation (Interior nodes)
		for (i = 1; i < Nx - 1; i++) {
			for (j = 1; j < Ny; j++) {
				u_E = 0.5 * (u(i, j) + u(i + 1, j));
				u_W = 0.5 * (u(i, j - 1) + u(i + 1, j - 1));
				v_N = 0.5 * (v(i - 1, j) + v(i, j));
				v_S = 0.5 * (v(i, j) + v(i + 1, j));

				a_E = -0.5 * u_E * h + nu;
				a_W = 0.5 * u_W * h + nu;
				a_N = -0.5 * v_N * h + nu;
				a_S = 0.5 * v_S * h + nu;

				a_n = 0.5 * u_E * h - 0.5 * u_W * h + 0.5 * v_N * h - 0.5 * v_S * h + 4 * nu;
				A_n = -h;

				d_n(i, j) = A_n / a_n;

				v_star(i, j) = (a_E * v(i, j + 1) + a_W * v(i, j - 1) + a_N * v(i - 1, j) + a_S * v(i + 1, j)) / a_n + d_n(i, j) * (p(i, j + 1) - p(i, j));
			}
		}

		// y-momentum BC's
		for (i = 0; i < Ny; i++) {
			v_star(i, 0) = -v_star(i, 1);
			v_star(i, Ny) = -v_star(i, Ny - 1);
		}

		for (j = 1; j < Nx; j++) {
			v_star(0, j) = 0;
			v_star(Nx - 1, j) = 0;
		}

		// Reset pressure corrections
		for (i = 0; i < Nx; i++) {
			for (j = 0; j < Ny; j++) {
				pc(i, j) = 0;
			}
		}

		// Continuity equation for the interior nodes (pressure correction)
		for (i = 1; i < Nx; i++) {
			for (j = 1; j < Ny; j++) {
				a_E = -d_e(i, j) * h;
				a_W = -d_e(i, j - 1) * h;
				a_N = -d_n(i - 1, j) * h;
				a_S = -d_n(i, j) * h;
				
				a_P = a_E + a_W + a_N + a_S;

				b(i, j) = -(u_star(i, j) - u_star(i, j - 1)) * h + (v_star(i, j) - v_star(i - 1, j)) * h;

				pc(i, j) = (a_E * pc(i, j + 1) + a_W * pc(i, j - 1) + a_N * pc(i - 1, j) + a_S * pc(i + 1, j) + b(i, j)) / a_P;
			}
		}

		// Correct the pressure field
		for (i = 1; i < Nx; i++) {
			for (j = 1; j < Ny; j++) {
				p_new(i, j) = p(i,j) + alpha_p*pc(i,j);
			}
		}

		// Continuity equation at the boundary
		for (j = 0; j < Ny; j++) {
			p_new(0, j) = p_new(1, j);
			p_new(Nx, j) = p_new(Nx - 1, j);
		}

		for (i = 1; i < Nx; i++) {
			p_new(i, 0) = p_new(i, 1);
			p_new(i, Ny) = p_new(i, Ny - 1);
		}

		// Correct u velocities
		for (i = 1; i < Nx; i++) {
			for (j = 1; j < Ny - 1; j++) {
				u_new(i, j) = u_star(i, j) + alpha * d_e(i, j) * (pc(i, j + 1) - pc(i, j));
			}
		}

		// Update x-momentum BC's
		for (j = 0; j < Ny; j++) {
			u_new(0, j) = 2 - u_new(1, j);
			u_new(Nx, j) = -u_new(Nx - 1, j);
		}

		for (i = 1; i < Nx; i++) {
			u_new(i, 0) = 0;
			u_new(i, Ny - 1) = 0;
		}

		// Correct v velocities
		for (i = 1; i < Nx - 1; i++) {
			for (j = 1; j < Ny; j++) {
				v_new(i, j) = v_star(i, j) + alpha * d_n(i, j) * (pc(i, j) - pc(i + 1, j));
			}
		}

		// Update y-momentum BC's
		for (i = 0; i < Nx; i++) {
			v_new(i, 0) = -v_new(i, 1);
			v_new(i, Ny) = -v_new(j, Nx - 1);
		}

		for (j = 1; j < Nx; j++) {
			v_new(0, j) = 0;
			v_new(Nx - 1, j) = 0;
		}

		// Continuity residual as error measure
		error = 0;
		for (i = 1; i < Nx; i++) {
			for (j = 1; j < Ny; j++) {
				error = error + pow(pow(b(i, j), 2), 0.5);
			}
		}

		// Update parameters for next itteration
		error = 0;
		for (i = 1; i < Nx; i++) {
			for (j = 1; j < Ny; j++) {
				u(i,j) = u_new(i,j);
				v(i,j) = v_new(i,j);
				p(i,j) = p_new(i,j);
			}
		}
		iterations = iterations + 1;
	}
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//unsigned long itime = 0;
	//unsigned long k = 0;
	//char fileName[50] = "solution_0.vtk";

	//initializeVel();
	//initializePres();

	//storeVTKSolution(u,v,p,G,fileName);
	///*========================================================================================
	// * assemble the LHS for the X momentum equation Au (already defined in the sturcture) 
	// *========================================================================================
 //	 */ 
	//assembleA('U');
	///*========================================================================================
	// * assemble the LHS for the Y momentum equation Av (already defined in the sturcture)
	// *========================================================================================
 //	 */
	//assembleA('V');

	//while(itime < Params.nTimeSteps){
	//	printf("========================================================\n");
	//	printf("Time step = %lu, Current Solution time = %f\n",itime+1, (itime+1)*Params.dt);
	//	printf("========================================================\n");
	//	// Fcu_ flux of convection term for u at n-1
	//	Fcu_ = Fcu;
	//    /*========================================================================================
	//     * given u and v, calculate udu/dx+vdu/dy in conservative form
	//     *========================================================================================
 //	     */
	//	computeConvectiveXFlux(Fcu, u, v, G);
	//    /*========================================================================================
	//     * given u and v, calculate d^2u/dx^2 with central difference scheme
	//     *========================================================================================
 //	     */
	//	computeDiffusiveFlux(Fdu, u, G);
	//    /*========================================================================================
	//     * assemble the residual without the pressure term for X momentum
	//     *========================================================================================
 //	     */
	//	assembleb0('U');
	//	// Fcv_ flux of convection term for v at n-1
	//	Fcv_ = Fcv;
	//    /*========================================================================================
	//     * given u and v, calculate udv/dx+vdv/dy in conservative form
	//     *========================================================================================
 //	     */
	//	computeConvectiveYFlux(Fcv, u, v, G);
	//    /*========================================================================================
	//     * given u and v, calculate d^2u/dx^2 with central difference scheme
	//     *========================================================================================
 //	     */
	//	computeDiffusiveFlux(Fdv, v, G);
	//    /*========================================================================================
	//     * assemble the residual without the pressure term for Y momentum
	//     *========================================================================================
 //	     */
	//	assembleb0('V');

	//	k = 0;

	//	//Set intermediate flow variables prior to outer iterations
	//	// Initialize predictors
	//	u_ = u; v_ = v; p_ = p;

	//	//Run outer iterations
	//	while(k < Params.maxOuterIterations){
	//		printf("----------------------Iteration = %lu-------------------------\n",k);

	//	/*========================================================================================
	//     * given p calculate pressure gradient
	//     *========================================================================================
 //	     */
	//		computeGradient(gradPx, gradPy, p_, G);
	//	/*========================================================================================
	//     * predict intermediate u_ and v_ according to guessed p_
	//     *========================================================================================
 //	     */			
	//		solveXMomPredictor();
	//		solveYMomPredictor();
	//	/*========================================================================================
	//     * solve the Poisson's equation for the pressure correction pc
	//     * For simplicity, solve it explicitly. No Newton-Raphson iteration is needed
	//     * - (1/a)\Delta p = - (du/dx+du/dy+d^4p/dx^4+d^4p/dy^4)
	//     * Note that we have handeled the boundary condition for you.
	//     * You just need to take care of internal nodes
	//     *========================================================================================
 //	     */
	//		solvePresCorrector();
	//	/*========================================================================================
	//     * update the pressure p_ according to pc
	//     * update the velocity u_,v_ 
	//     * (note that u_ and v_ will be overwritten in the next outer iteration.
	//     * u_ and v_ will be passed to the next step when the outer iteration converges.)
	//     *========================================================================================
 //	     */
	//		correctFlowVars();
	//		k++;
	//		printf("--------------------------------------------------------------\n");
	//		if(checkConvergence())
	//			break;
	//	}

	//	if(k < Params.maxOuterIterations)
	//		printf("pressure-velocity coupling converged in %lu iterations.\n",k);
	//	else
	//		printf("pressure-velocity coupling limited to %lu iterations.\n",k);

	//	if((itime+1) % 50 == 0){
	//		sprintf(fileName,"solution_%lu.vtk",itime+1);
	//		storeVTKSolution(u,v,p,G,fileName);
	//	}

	//	//Update velocity and pressure for next time step
	//	u = u_; v = v_; p = p_;
	//	itime++;
	//}

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

void NSDriver::assembleb0(char flowVar)
{
	switch(flowVar)
	{
		case 'U':   // Added explicit block {}
		{
			int bu0 = 0;   // Added int
			break;
		}
		case 'V':   // Added explicit block {}
		{
			int bv0 = 0;   // Added int
			break;
		}
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
			}	
		break;

		case 'V':
		for(size_t i = 1; i < G.Nx()-1; i++)
			for(size_t j = 1; j < G.Ny()-1; j++)
			{
				Av(i,j,0) = 0;
			}	
		break;
	}
}


void NSDriver::assembleAp()
{
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
		
			Ap(1,j,0) = 0;
			Ap(1,j,1) = 0;
			Ap(1,j,2) = 0;
			Ap(1,j,3) = 0;
			Ap(1,j,4) = 0;
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
			bp(i,j) = 0 ;

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