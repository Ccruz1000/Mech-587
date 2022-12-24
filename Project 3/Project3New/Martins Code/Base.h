#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstdlib>

#ifndef BASE_H
#define BASE_H
#define INITVEC true
#define TOL 1e-8
#define PI 4.0*atan(1.0)
#define MIN(A,B) ((A) < (B)? (A):(B))
#define MAX(A,B) ((A) > (B)? (A):(B))

/**
 *	The Operations class contains the underlying operations such as memory allocation, and deallocation. 
 *  It also contains vector operations such as copying, addition, scalar multiplication, dot products, and Norms of vectors.
 */


class Operations{
public:
	void alloc1D(double **V, const unsigned long &N, bool isInit = true, double def = 0.0) const;
	void deAlloc1D(double **V, const unsigned long &N) const;
	void copyArray(double *S, double *T, const unsigned long &N) const;
	void addVec(double *R, const double *V1, const double *V2, const unsigned long &N) const;
	void subtractVec(double *R, const double *V1, const double *V2, const unsigned long &N) const;
	void scaleVec(double *R, const double &a, const double *V1, const unsigned long &N) const;
	void dotVec(double &R, const double *V1, const double *V2, const unsigned long &N) const;
	void findAbsMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const;
	void findMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const;
	void findMin(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const;
};



class Grid;
class Vector;
class Matrix;
class Config;
class Solver;	//Solves the Linear System Ax = b


/**
 * The Grid class initializes the 2-Dimensional Cartesian Grid, necessary for handling the 
 * problem. 
 * 
 * It contains the following data members:
 * X, Y 	-> Co-ordinates of Grid in X and Y directions
 * nx, ny 	-> No. of grid-points in X and Y directions.
 * dX, dY	-> Grid spacing in X and Y directions
 * xlim[2], ylim[2] -> Domain Limits in X and Y directions
 * N 		-> Total no. of grid-points in the Grid, which will be the same as the size of the Solution vector.
 * 
 * 
 * The methods for the Grid class are:
 * Grid()			-> Default constructor
 * Grid(Nx, Ny)  	-> Constructor which takes in no. of Grid Points as input, and initializes a unit square as default domain.
 * Grid(Nx, Ny, xlim, ylim) -> Constructor which takes in no. of Grid Points and domain limits as input 
 * x(i)			-> Returns the x-co-ordinate of ith point along x-direction
 * y(j)			-> Returns the y-co-ordinate of jth point along y-direction
 *  
 */

class Grid
{
	private:
	double *X, *Y;
	unsigned long nx, ny;
	double dX, dY;
	unsigned long N;
	double xlim[2], ylim[2];
	Operations O;
	
	public:
	Grid();
	Grid(unsigned long Nx, unsigned long Ny);
	Grid(unsigned long Nx, unsigned long Ny, double xlim[2], double ylim[2]);
	~Grid();
	void setGridSize(const unsigned long Nx, const unsigned long Ny);
	inline double x(const unsigned long i) const {return X[i];}
	inline double y(const unsigned long j) const {return Y[j];}
	inline double dx() const {return dX;}
	inline double dy() const {return dY;}
	inline size_t Nx() const {return nx;}
	inline size_t Ny() const {return ny;}
	inline size_t size() const {return N;}
	void storeVTKStructured(const Vector &u, const char *fname) const;
};

class Matrix
{
	unsigned long N, Nx, Ny;
	double *A[5];
	bool isInit;
	Operations O;

	bool isValid(unsigned short pos, unsigned long idxNode) const;
public:
	Matrix();
	Matrix(unsigned long iNx, unsigned long iNy);
	~Matrix();
	inline size_t size() const {return N;}
	inline size_t sizeNx() const {return Nx;}
	inline size_t sizeNy() const {return Ny;}
	void setSize(unsigned long iNx, unsigned long iNy);
	double& operator()(unsigned long i, unsigned long j, unsigned short pos);
	double operator()(unsigned long i, unsigned long j, unsigned short pos) const;

	double& operator()(unsigned long i, unsigned short pos);
	double operator()(unsigned long i, unsigned short pos) const;

	void storeA(char filename[50]);


};

class Vector{
	unsigned long N, Nx, Ny;
	double *b;
	bool isInit;
	Operations O;

	inline bool isValid(unsigned long idxNode) const {return idxNode < N ? true:false; }
public:
	Vector();
	Vector(unsigned long Nx, unsigned long Ny, bool isInit = true, double initVal = 0.0);
	Vector(const Vector &V1);
	~Vector();
	inline size_t size() const{return N;}
	inline size_t sizeNx() const{return Nx;}
	inline size_t sizeNy() const{return Ny;}
	void setSize(unsigned long Nx, unsigned long Ny);
	double& operator()(unsigned long i, unsigned long j);
	double operator()(unsigned long i, unsigned long j) const;

	double& operator()(unsigned long i);
	double operator()(unsigned long i) const;

	Vector operator -();
	Vector operator + (const Vector &V1) const;
	Vector operator - (const Vector &V1) const;
	Vector operator * (const double &S) const;
	Vector operator + (const Vector &V1);
	Vector operator - (const Vector &V1);
	Vector operator * (const double &S);
	Vector operator = (const Vector &V1);
	friend const double operator % (const Vector &V1, const Vector &V2);	//For dot product
	friend Vector operator * (double const &S, Vector const &V);
	friend Vector operator * (const Matrix &A, Vector const &v);

	const double L2Norm();
	const double LinfNorm(unsigned long &i, unsigned long &j);
	const double GetMax();
	const double GetMin();

	void storeV(char filename[50]);
};

class Solver
{
	private:
		const unsigned short solverCode;
		const unsigned short solverPreConditioner;
		const unsigned long maxIterations;
		const double omegaSOR;
		const double alphaS;
		const double relResidual;
		const double absResidual;

		size_t Nx, Ny;
		Vector R0, R, v, p, h, s, t, y, w, z;
		double rhoi, rho0, alpha, omega, beta;
		Matrix LU;
		

	void solveGS(Vector &u, const Matrix &A, const Vector &b);
	void decomposeILU(const Matrix &A);
	void solveLU(Vector &u, const Matrix &LU, const Vector &b);
	void solveBiCG(Vector &u, const Matrix &A, const Vector &b);
	
	public:
	Solver();
	~Solver();
	void solve(Vector&u, const Matrix &A, const Vector&b);
	
};

// Add a configuration file reader class, and a configuration class
// The configuration class will store the parameters for the problem
// Reynolds Number, Grid Size, Linear Solver type, Linear Solver Preconditioner
// The Config class will be a global singleton. Only one such object can be brought into existence,
// and will be shared across the solver. Further, it cannot be modified over the course of the program.

class Config
{
	// Grid Parameters
	unsigned long nx, ny;		//Grid points along x and y directions

	// Navier-Stokes Parameters
	double Dt;					//Time-step size.
	double re;					//Reynolds Number.
	unsigned long maxTimeSteps;	//Max-time steps
	double alphaP;				//Pressure under-relaxation factor
	double alphaU;				//Momentum under-relaxation factor
	double convergence;			//Convergence criteria for NS equations
	unsigned long maxOuterIterations;	//Maximum number of outer iterations

	//Linear Solver Parameters
	unsigned short linSolver;		// Choice of linear Solver - 1 -> Gauss Seidel, 2 ->  BiCGStab
	double omega;					// Over-relaxation parameter
	double relResidual;				// Convergence criteria for relative residual in Linear Solver
	double absResidual;				// Convergence criteria for absolute residual in Linear Solver
	unsigned long maxIterations;	// Max Iterations of Linear Solver
	
	// Choice of pre-conditioner - Only for BiCGStab. 0 -> No preconditioning, 1 -> ILU, 2 -> Stone's Preconditioner
	unsigned short preConditioner;	
	double alphaS;					// Stone's parameter for preconditioning


	//File reader
	FILE *inputfile;

	//Private Constructor 
	Config();							// Private Constructor
	~Config(){}						// Private distructor
	Config(const Config &C);			// Cannot create a copy
	Config& operator=(const Config&);	// Cannot assign

	static Config *instance;
public:
	static Config *Instance()
	{
		if(!instance)
			instance = new Config();
		return instance;
	}

	inline const size_t Nx() const {return nx;}
	inline const size_t Ny() const {return ny;}

	inline const double dt() const {return Dt;}
	inline const double Re() const {return re;}
	inline const unsigned long nTimeSteps() const {return maxTimeSteps;}
	inline const double urfPres() const {return alphaP;}
	inline const double urfMom() const {return alphaU;}
	inline const double conPresMom() const {return convergence;}
	inline const unsigned short maxOuterIter() const {return maxOuterIterations;}

	inline const unsigned short selectlinSolver() const {return linSolver;}
	inline const double orfSOR() const {return omega;}
	inline const double solverRelResidual() const {return relResidual;}
	inline const double solverAbsResidual() const {return absResidual;}
	inline const unsigned long maxSolverIterations() const{return maxIterations;}
	inline const unsigned short selectPreConditioner() const {return preConditioner;}
	inline const double solverPreCondFactor() const {return alphaS;}
};

void storeVTKStructured(const Vector &u, const Grid &G, const char *fileName);
void storeVTKSolution(const Vector &u, const Vector &v, const Vector &p, const Grid &G, const char* fileName) ;


#endif