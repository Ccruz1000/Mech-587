#include "Base.h"
#include "Convection_Diffusion.h"

#ifndef NSDRIVER_H
#define NSDRIVER_H

struct NSParameters
{
	const double dt;			//Time-step size
	const double Re;			//Reynolds Number.
	const unsigned long nTimeSteps;	//Max-time steps
	const double alphaP;				//Pressure under-relaxation factor
	const double alphaU;				//Momentum under-relaxation factor
	const double convergence;			//Convergence criteria for NS equations
	const unsigned long maxOuterIterations;	//Maximum number of outer iterations

	NSParameters();
	~NSParameters(){}
};


class NSDriver
{
	const Grid G;
	const NSParameters Params;
	Solver S;

	//Solution vectors for velocity and pressure
	Vector u, v, p;	 

	//Vectors for velocity predictor
	Vector u_, v_, p_;
	Vector gradPx, gradPy;
	Vector uf, vf;	//Face velocities

	//Corrector quantities
	Vector pc;
	Vector momResx, momResy;

	//Matrices for LHS 
	Matrix Au, Av, Ap;	

	//R.H.S vectors for momentum and pressure equations
	Vector bu, bv, bp;
	Vector bu0, bv0;

	//Flux vectors for convection and Diffusion
	Vector Fcu, Fcv, Fcu_, Fcv_;
	Vector Fdu, Fdv;

	void assembleA(char flowVar);
	// void assembleb0(char flowVar);
	void assembleb0(char flowVar,const Vector &Fc, const Vector &Fc_,const Vector &Fd, const Vector &vel, unsigned long itime);
	void applyMomBC(char flowVar);

	void assembleFaceVel();
	void assembleAp();
	void assemblebp();
	void anchorpc();
	void applyPresBC();

	void initializeVel();
	void initializePres();

	void solveXMomPredictor();
	void solveYMomPredictor();
	void enforceZeroDivergence();
	void solvePresCorrector();
	void correctFlowVars();
	bool checkConvergence();

public:
	NSDriver();
	~NSDriver();
	void Run();
	void TestFlux();

};

#endif


