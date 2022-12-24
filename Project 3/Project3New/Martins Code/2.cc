#include<cmath>
#include "Convection_Diffusion.h"

// Setup the functions to solve the convection diffusion equations. 
void computeConvectiveXFlux(Vector &Fc, const Vector &u, const Vector &v, const Grid &G)
{

	unsigned long i,j;
	double dx, dy;
	dx = G.dx(); dy = G.dy();

	double ue, uw, un, us;
	double ve, vw, vn, vs;

	ue = 0.5 * (u(i + 1, j) + u(i, j));
	uw = 0.5 * (u(i - 1, j) + u(i, j));
	un = 0.5 * (u(i, j + 1) + u(i, j));
	us = 0.5 * (u(i, j - 1) + u(i, j));

	ve = 0.5 * (v(i + 1, j) + v(i, j));
	vw = 0.5 * (v(i - 1, j) + v(i, j));
	vn = 0.5 * (v(i, j + 1) + v(i, j));
	vs = 0.5 * (v(i, j - 1) + v(i, j));


	for(size_t i = 1; i < G.Nx()-1; i++)
		for(size_t j = 1; j < G.Ny()-1; j++){

			Fc(i, j) = (((ue * ue) - (uw * uw)) / dx) + (((un * vn) - (us * vs)) / dy);
		}
}

void computeConvectiveYFlux(Vector &Fc, const Vector &u, const Vector &v, const Grid &G)
{

	unsigned long i,j;
	double dx, dy;
	dx = G.dx(); dy = G.dy();

	double ue, uw, un, us;
	double ve, vw, vn, vs;

	ue = 0.5 * (u(i + 1, j) + u(i, j));
	uw = 0.5 * (u(i - 1, j) + u(i, j));
	un = 0.5 * (u(i, j + 1) + u(i, j));
	us = 0.5 * (u(i, j - 1) + u(i, j));

	ve = 0.5 * (v(i + 1, j) + v(i, j));
	vw = 0.5 * (v(i - 1, j) + v(i, j));
	vn = 0.5 * (v(i, j + 1) + v(i, j));
	vs = 0.5 * (v(i, j - 1) + v(i, j));

	for(size_t i = 1; i < G.Nx()-1; i++)
		for(size_t j = 1; j < G.Ny()-1; j++){

			Fc(i, j) = (((vn * vn) - (vs * vs)) / dy) + (((ue * ve) - (uw * vw)) / dx);
		}
}

void computeDiffusiveFlux(Vector &Fd, const Vector &phi, const Grid &G)
{
	unsigned long i,j;
	double dx, dy;
	dx = G.dx(), dy = G.dy();

	double dphidxe, dphidxw, dphidyn, dphidys;

	for(size_t i = 1; i < G.Nx()-1; i++)
		for(size_t j = 1; j < G.Ny()-1; j++){

			Fd(i, j) = (0.01) * (((phi(i - 1, j) - 2 * phi(i, j) + phi(i + 1, j)) / (dx * dx)) + ((phi(i, j - 1) - 2 * phi(i, j) + phi(i, j + 1)) / (dy * dy)));
		}
}

void computeGradient(Vector &gradphix, Vector &gradphiy, const Vector &phi, const Grid &G)
{
	unsigned long i,j;
	double 			dx;			//x- grid point spacing
	double 			dy;			//y- grid point spacing
	dx = G.dx(); dy = G.dy();

	double pe, pw, pn, ps;

	pe = 0.5 * ((phi(i + 1, j)) + (phi(i, j)));
	pw = 0.5 * ((phi(i - 1, j)) + (phi(i, j)));
	pn = 0.5 * ((phi(i, j + 1)) + (phi(i, j)));
	ps = 0.5 * ((phi(i, j - 1)) + (phi(i, j)));

	//Construct gradient using CDS over interior points
	for(size_t i = 1; i < G.Nx() -1; i++)
		for(size_t j = 1; j < G.Ny()-1; j++)
		{
			gradphix(i, j) = (pe - pw) / dx;
			gradphiy(i, j) = (pn - ps) / dy;
		}
}


// void testFlux(const Grid &G)
// {
// 	unsigned long Nx, Ny;
// 	Nx = G.Nx(), Ny = G.Ny();
// 	double dy = G.dy();
// 	double dx = G.dx();

// 	double x[4], y[4];
// 	double u[4], v[4], p[4];
// 	double dudx[4], dudy[4];
// 	double dvdx[4], dvdy[4];
// 	double dpdx[4], dpdy[4];

// 	Vector ue(Nx,Ny), ve(Nx,Ny), pe(Nx,Ny);
// 	Vector FcuExact(Nx,Ny), FduExact(Nx,Ny), FpuExact(Nx,Ny);
// 	Vector FcvExact(Nx,Ny), FdvExact(Nx,Ny), FpvExact(Nx,Ny);

// 	Vector Fcu(Nx,Ny), Fcv(Nx,Ny), Fdu(Nx,Ny), Fdv(Nx,Ny), Fpu(Nx,Ny), Fpv(Nx,Ny);

// 	for(size_t i = 1; i < Nx-1; i++)
// 		for(size_t j = 1; j < Ny-1; j++){
// 			x[0] = G.x(i) + 0.5*dx; y[0] = G.y(j);
// 			x[1] = G.x(i); y[1] = G.y(j) + 0.5*dy;
// 			x[2] = G.x(i) - 0.5*dx; y[2] = G.y(j);
// 			x[3] = G.x(i); y[3] = G.y(j) - 0.5*dy;

// 			ue(i,j) = 1.0*sin(PI*G.x(i))*sin(2*PI*G.y(j));
// 			ve(i,j) = 1.0*sin(2*PI*G.x(i))*sin(PI*G.y(j));
// 			pe(i,j) = 1.0*cos(PI*G.x(i))*cos(PI*G.y(j));

// 			for(int k = 0; k < 4; k++)
// 			{
// 				u[k] = 1.0*sin(PI*x[k])*sin(2*PI*y[k]);
// 				v[k] = 1.0*sin(2*PI*x[k])*sin(PI*y[k]);
// 				p[k] = 1.0*cos(PI*x[k])*cos(PI*y[k]);

// 				dudx[k] = 1.0*PI*cos(PI*x[k])*sin(2*PI*y[k]);
// 				dudy[k] = 1.0*2*PI*sin(PI*x[k])*cos(2*PI*y[k]);

// 				dvdx[k] = 1.0*2*PI*cos(2*PI*x[k])*sin(PI*y[k]);
// 				dvdy[k] = 1.0*PI*sin(2*PI*x[k])*cos(PI*y[k]);

// 				dpdx[k] = -1.0*PI*sin(PI*x[k])*cos(PI*y[k]);
// 				dpdy[k] = -1.0*PI*cos(PI*x[k])*sin(PI*y[k]);
// 			}

// 			FcuExact(i,j) = (u[0]*u[0] - u[2]*u[2])*dy + (v[1]*u[1] - v[3]*u[3])*dx;
// 			FduExact(i,j) = (dudx[0] - dudx[2])*dy + (dudy[1] - dudy[3])*dx;
// 			FpuExact(i,j) = (p[0] - p[2])*dy;

// 			FcvExact(i,j) = (u[0]*v[0] - u[2]*v[2])*dy + (v[1]*v[1] - v[3]*v[3])*dx;
// 			FdvExact(i,j) = (dvdx[0] - dvdx[2])*dy + (dvdy[1] - dvdy[3])*dx;
// 			FpvExact(i,j) = (p[1] - p[3])*dy;			
// 		}

// 	for(size_t i = 0; i < Nx; i++)
// 		for(size_t j = 0; j<=Ny-1; j+=Ny-1){
// 			ue(i,j) = 1.0*sin(PI*G.x(i))*sin(2*PI*G.y(j));
// 			ve(i,j) = 1.0*sin(2*PI*G.x(i))*sin(PI*G.y(j));
// 			pe(i,j) = 1.0*cos(PI*G.x(i))*cos(PI*G.y(j));	
// 		}

// 	for(size_t j = 0; j < Ny; j++)
// 		for(size_t i = 0; i<=Nx-1; i+=Nx-1){
// 			ue(i,j) = 1.0*sin(PI*G.x(i))*sin(2*PI*G.y(j));
// 			ve(i,j) = 1.0*sin(2*PI*G.x(i))*sin(PI*G.y(j));
// 			pe(i,j) = 1.0*cos(PI*G.x(i))*cos(PI*G.y(j));	
// 		}

// 	//Compute convective, diffusive flux and pressure gradient
// 	computeConvectiveFlux(Fcu,ue,ve,ue,G);
// 	computeConvectiveFlux(Fcv,ue,ve,ve,G);
// 	computeDiffusiveFlux(Fdu,ue,G);
// 	computeDiffusiveFlux(Fdv,ve,G);
// 	computeGradient(Fpu,Fpv,pe,G);
	
// 	FcuExact = 0.5*FduExact - FcuExact;
// 	FcvExact = 0.5*FdvExact - FcvExact;

// 	Fcu = 0.5*Fdu - Fcu;
// 	Fcv = 0.5*Fdv - Fcv;

// 	FcuExact = FcuExact - Fcu;FcvExact = FcvExact - Fcv;
// 	printf("Convection_Diffusion Test executed for x-Mom. Flux Error = %14.12e %lu\n",FcuExact.L2Norm(),FcuExact.size());
// 	printf("Convection_Diffusion Test executed for y-Mom. Flux Error = %14.12e\n",FcvExact.L2Norm());

// 	FpuExact = FpuExact - Fpu;
// 	FpvExact = FpvExact - Fpv;
// 	printf("Pressure Test executed for x-Mom. Flux Error = %14.12e\n",FpuExact.L2Norm());
// 	printf("Pressure Test executed for y-Mom. Flux Error = %14.12e\n",FpvExact.L2Norm());

// }

// void Driver::testDiffusiveFlux()
// {
// 	unsigned long Nx, Ny;
// 	Nx = G->Nx(), Ny = G->Ny();
// 	double dx, dy;
// 	dx = G->dx(), dy = G->dy();

// 	double Re = C->Re();
// 	double due, duw, dun, dus;
// 	double dve, dvw, dvn, dvs;

// 	Vector fDExact(Nx,Ny);
// 	computeDiffusiveFlux(FlowVar::U);
// 	for(unsigned long i = 1; i < Nx-1; i++)
// 		for(unsigned long j = 1; j < Ny-1; j++){
// 			due = 0.0, duw = 0.0;
// 			dun = (1.0/Re)*(G->y(j+1) - G->y(j))/dy;
// 			dus = (1.0/Re)*(G->y(j) - G->y(j-1))/dy;

// 			fDExact(i,j) = (due-duw)*dy + (dun-dus)*dx;
// 		}
// 	fDExact = fDExact-fV.uDn;
// 	printf("Diffusion Test executed for x-Mom. Flux Error = %14.12e\n",fDExact.L2Norm());

// 	computeDiffusiveFlux(FlowVar::V);
// 	for(unsigned long i = 1; i < Nx-1; i++)
// 		for(unsigned long j = 1; j < Ny-1; j++){
// 			dvn = 0.0, dvs = 0.0;
// 			dve = (1.0/Re)*(G->x(i) - G->x(i+1))/dx;
// 			dvw = (1.0/Re)*(G->x(i-1) - G->x(i))/dx;

// 			fDExact(i,j) = (dve-dvw)*dy + (dvn-dvs)*dx;
// 		}
// 	fDExact = fDExact-fV.vDn;
// 	printf("Diffusion Test executed for y-Mom. Flux Error = %14.12e\n",fDExact.L2Norm());
// }

// void Driver::testPressureGradient()
// {
// 	unsigned long Nx, Ny;
// 	Nx = G->Nx(), Ny = G->Ny();
// 	double dx, dy;
// 	dx = G->dx(), dy = G->dy();

// 	double dpdx, dpdy;

// 	Vector gradPxE(Nx,Ny);
// 	Vector gradPyE(Nx,Ny);
// 	computePressureGradient(FlowVar::U);
// 	computePressureGradient(FlowVar::V);
// 	for(unsigned long i = 1; i < Nx-1; i++)
// 		for(unsigned long j = 1; j < Ny-1; j++){
// 			gradPxE(i,j) = -(pow((0.5*(G->x(i+1) + G->x(i))),1) - pow((0.5*(G->x(i) + G->x(i-1))),1))*dy;
// 			gradPyE(i,j) = -(pow((0.5*(G->y(j+1) + G->y(j))),1) - pow((0.5*(G->y(j) + G->y(j-1))),1))*dx;
// 		}
// 	gradPxE = gradPxE - fV.gradPx;
// 	gradPyE = gradPyE - fV.gradPy;

// 	printf("Pressure Test executed. gradPx error = %14.12e, gradPy error = %14.12e\n",gradPxE.L2Norm(), gradPyE.L2Norm());



// Driver::Driver(const Grid *iG, const Config *iC) : G(iG), C(iC),
// 					Au(iG->Nx(),iG->Ny()), bu(iG->Nx(),iG->Ny()),
// 					Av(iG->Nx(),iG->Ny()), bv(iG->Nx(),iG->Ny()),
// 					Ap(iG->Nx(),iG->Ny()), bp(iG->Nx(),iG->Ny()),
// 					u_(iG->Nx(),iG->Ny()), v_(iG->Nx(),iG->Ny()),
// 					p_(iG->Nx(),iG->Ny())
// {
// 	initializeFlowVars(fV, *G);
// }

// Driver::~Driver()
// {
// 	G = NULL;
// 	C = NULL;
// }

// void Driver::computeConvectiveFlux(unsigned short flowEq)
// {
// 	//Localize Variables
// 	double 			Vc;			//Cell volume
// 	double 			dx;			//x- grid point spacing
// 	double 			dy;			//y- grid point spacing

// 	double			ue, uw, un, us;
// 	double 			ve, vw, vn, vs;
// 	double			fcx, fcy;

// 	dx = G->dx(); dy = G->dy();

// 	switch(flowEq)
// 	{
// 		case FlowVar::U:
// 			//Copy convective flux at time level n to time level n-1
// 			fV.uCn_ = fV.uCn;

// 			//Construct Convective flux using CDS over interior points
// 			for(long i = 1; i < G->Nx()-1; i++)
// 				for(long j = 1; j < G->Ny()-1; j++)
// 				{
// 					ue = 0.5*(fV.u(i+1,j) + fV.u(i,j));
// 					uw = 0.5*(fV.u(i,j) + fV.u(i-1,j));

// 					un = 0.5*(fV.u(i,j+1) + fV.u(i,j));
// 					us = 0.5*(fV.u(i,j) + fV.u(i,j-1));
// 					vn = 0.5*(fV.v(i,j+1) + fV.v(i,j));
// 					vs = 0.5*(fV.v(i,j) + fV.v(i,j-1));

// 					fcx = (ue*ue - uw*uw)*dy;
// 					fcy = (un*vn - us*vs)*dx;

// 					fV.uCn(i,j) = fcx + fcy;
// 				}
// 			break;

// 		case FlowVar::V:
// 			//Copy convective flux at time level n to time level n-1
// 			fV.vCn_ = fV.vCn;

// 			//Construct Convective flux using CDS over interior points
// 			for(long i = 1; i < G->Nx()-1; i++)
// 				for(long j = 1; j < G->Ny()-1; j++)
// 				{
// 					ue = 0.5*(fV.u(i+1,j) + fV.u(i,j));
// 					uw = 0.5*(fV.u(i,j) + fV.u(i-1,j));
// 					ve = 0.5*(fV.v(i+1,j) + fV.v(i,j));
// 					vw = 0.5*(fV.v(i,j) + fV.v(i-1,j));
					
// 					vn = 0.5*(fV.v(i,j+1) + fV.v(i,j));
// 					vs = 0.5*(fV.v(i,j) + fV.v(i,j-1));

// 					fcx = (ue*ve - uw*vw)*dy;
// 					fcy = (vn*vn - vs*vs)*dx;

// 					fV.vCn(i,j) = fcx + fcy;
// 				}
// 			break;
// 	}
// }

// void Driver::computeDiffusiveFlux(unsigned short flowEq)
// {
// 	double 			Vc;			//Cell volume
// 	double 			dx;			//x- grid point spacing
// 	double 			dy;			//y- grid point spacing

// 	double			due, duw, dun, dus;
// 	double 			dve, dvw, dvn, dvs;
// 	double			fdx, fdy;

// 	double			Re;

// 	dx = G->dx(); dy = G->dy();
// 	Re = C->Re();

// 	switch(flowEq)
// 	{
// 		case FlowVar::U:
// 			//Construct Diffusive flux using CDS over interior points
// 			for(long i = 1; i < G->Nx()-1; i++)
// 				for(long j = 1; j < G->Ny()-1; j++)
// 				{
// 					due = (fV.u(i+1,j) - fV.u(i,j))/dx;
// 					duw = (fV.u(i,j) - fV.u(i-1,j))/dx;
// 					dun = (fV.u(i,j+1) - fV.u(i,j))/dy;
// 					dus = (fV.u(i,j) - fV.u(i,j-1))/dy;

// 					fdx = (due - duw)*dy;
// 					fdy = (dun - dus)*dx;

// 					fV.uDn(i,j) = (fdx + fdy)/Re;
// 				}
// 			break;

// 		case FlowVar::V:
// 			//Construct Diffusive flux using CDS over interior points
// 			for(long i = 1; i < G->Nx()-1; i++)
// 				for(long j = 1; j < G->Ny()-1; j++)
// 				{
// 					dve = (fV.v(i+1,j) - fV.v(i,j))/dx;
// 					dvw = (fV.v(i,j) - fV.v(i-1,j))/dx;
// 					dvn = (fV.v(i,j+1) - fV.v(i,j))/dy;
// 					dvs = (fV.v(i,j) - fV.v(i,j-1))/dy;

// 					fdx = (dve - dvw)*dy;
// 					fdy = (dvn - dvs)*dx;

// 					fV.vDn(i,j) = (fdx + fdy)/Re;
// 				}
// 			break;
// 	}
// }



// void Driver::computeA(unsigned short flowEq)
// {
// 	double dx; 
// 	double dy;
// 	double Vc;

// 	double dt, Re;

// 	double ap, ae, aw, an, as;

// 	dx = G->dx(); dy = G->dy();
// 	Vc = dx*dy;
// 	dt = C->dt();
// 	Re = C->Re();

// 	switch(flowEq)
// 	{
// 		case FlowVar::U:
// 			//Compute Transient term
// 			ap = 1.0*Vc/dt;
// 			//Compute Diffusion term
// 			ae = -0.5*(dy/dx)*(1.0/Re); aw = ae;
// 			an = -0.5*(dx/dy)*(1.0/Re); as = an;
// 			ap = ap - (ae + aw + an + as);
// 			//Save to Matrix
// 			for(long i = 1; i < G->Nx()-1; i++)
// 				for(long j = 1; j < G->Ny()-1; j++){
// 					Au(i,j,0) = aw;
// 					Au(i,j,1) = as;
// 					Au(i,j,2) = ap;
// 					Au(i,j,3) = an;
// 					Au(i,j,4) = ae;	
// 				}


// 		case FlowVar::V:
// 			//Compute Transient term
// 			ap = 1.0*Vc/dt;
// 			//Compute Diffusion term
// 			ae = -0.5*(dy/dx)*(1.0/Re); aw = ae;
// 			an = -0.5*(dx/dy)*(1.0/Re); as = an;
// 			ap = ap - (ae + aw + an + as);
// 			//Save to Matrix
// 			for(long i = 1; i < G->Nx()-1; i++)
// 				for(long j = 1; j < G->Ny()-1; j++){
// 					Av(i,j,0) = aw;
// 					Av(i,j,1) = as;
// 					Av(i,j,2) = ap;
// 					Av(i,j,3) = an;
// 					Av(i,j,4) = ae;	
// 				}
// 	}
// }

// void Driver::computeb0(unsigned short flowEq)
// {
// 	double Vc;
// 	double dx, dy;
// 	double dt;

// 	dx = G->dx(), dy = G->dy();
// 	Vc = dx*dy;
// 	dt = C->dt();

// 	switch(flowEq)
// 	{
// 		case FlowVar::U:
// 			//Add Transient Term
// 			bu0 = (Vc/dt)*fV.u + 0.5*fV.uDn;
// 			//Add ABS-2 term for Convection
// 			if(fV.timeId < 5)
// 				bu0 = bu0 - fV.uCn;
// 			else
// 				bu0 = bu0 - (1.5*fV.uCn - 0.5*fV.uCn_);
// 		case FlowVar::V:
// 			//Add Transient Term
// 			bv0 = (Vc/dt)*fV.v + 0.5*fV.vDn;
// 			//Add ABS-2 term for Convection
// 			if(fV.timeId < 5)
// 				bv0 = bv0 - fV.vCn;
// 			else
// 				bv0 = bv0 - (1.5*fV.vCn - 0.5*fV.vCn_);
// 	}
// }

// void Driver::setMomBC(unsigned short flowEq)
// {
// 	size_t Nx, Ny;
// 	Nx = G->Nx(); Ny = G->Ny();
// 	long i, j;

// 	//x-Momentum equation
// 	if(flowEq == FlowVar::U)
// 	{	
// 		j = 0;
// 		for(long i = 0; i < Nx; i++)
// 		{
// 			Au(i,j,0) = 0, Au(i,j,1) = 0, Au(i,j,2) = 1.0, Au(i,j,3) = 0.0, Au(i,j,4) = 0.0;
// 			u_(i,j) = 0.0; bu(i,j) = 0.0;
// 		}

// 		j = Ny-1;
// 		for(long i = 0; i < Nx; i++)
// 		{
// 			Au(i,j,0) = 0, Au(i,j,1) = 0, Au(i,j,2) = 1.0, Au(i,j,3) = 0.0, Au(i,j,4) = 0.0;
// 			u_(i,j) = 1.0; bu(i,j) = 1.0;
// 		}

// 		i = 0;
// 		for(long j = 0; j < Ny; j++)
// 		{
// 			Au(i,j,0) = 0, Au(i,j,1) = 0, Au(i,j,2) = 1.0, Au(i,j,3) = 0.0, Au(i,j,4) = 0.0;
// 			u_(i,j) = 0.0; bu(i,j) = 0.0;
// 		}

// 		i = Nx-1;
// 		for(long j = 0; j < Ny; j++)
// 		{
// 			Au(i,j,0) = 0, Au(i,j,1) = 0, Au(i,j,2) = 1.0, Au(i,j,3) = 0.0, Au(i,j,4) = 0.0;
// 			u_(i,j) = 0.0; bu(i,j) = 0.0;
// 		}
// 	}

// 	if(flowEq == FlowVar::V)
// 	{	
// 		j = 0;
// 		for(long i = 0; i < Nx; i++)
// 		{
// 			Av(i,j,0) = 0, Av(i,j,1) = 0, Av(i,j,2) = 1.0, Av(i,j,3) = 0.0, Av(i,j,4) = 0.0;
// 			v_(i,j) = 0.0; bv(i,j) = 0.0;
// 		}

// 		j = Ny-1;
// 		for(long i = 0; i < Nx; i++)
// 		{
// 			Av(i,j,0) = 0, Av(i,j,1) = 0, Av(i,j,2) = 1.0, Av(i,j,3) = 0.0, Av(i,j,4) = 0.0;
// 			v_(i,j) = 0.0; bv(i,j) = 0.0;
// 		}

// 		i = 0;
// 		for(long j = 0; j < Ny; j++)
// 		{
// 			Av(i,j,0) = 0, Av(i,j,1) = 0, Av(i,j,2) = 1.0, Av(i,j,3) = 0.0, Av(i,j,4) = 0.0;
// 			v_(i,j) = 0.0; bv(i,j) = 0.0;
// 		}

// 		i = Nx-1;
// 		for(long j = 0; j < Ny; j++)
// 		{
// 			Av(i,j,0) = 0, Av(i,j,1) = 0, Av(i,j,2) = 1.0, Av(i,j,3) = 0.0, Av(i,j,4) = 0.0;
// 			v_(i,j) = 0.0; bv(i,j) = 0.0;
// 		}
// 	}

// }

// void Driver::solveXMom()
// {
// 	//Compute pressure gradient
// 	computePressureGradient(FlowVar::U);

// 	//Compute b
// 	bu = bu0 + fV.gradPx;

// 	//Apply Boundary Conditions
// 	setMomBC(FlowVar::U);

// 	//Solve for u*
// 	solve(u_,Au,bu);
// }

// void Driver::solveYMom()
// {
// 	//Compute pressure gradient
// 	computePressureGradient(FlowVar::V);

// 	//Compute v
// 	bv = bv0 + fV.gradPy;

// 	//Apply Boundary Conditions
// 	setMomBC(FlowVar::V);

// 	//Solve for v*
// 	solve(v_,Av,bv);
// }


// void Driver::RunSolution()
// {
// 	int m = 0;
// 	//Since we are only working with diffusion terms implicitly, the matrices can
// 	//be assembled before the time stepping starts
// 	computeA(FlowVar::U);
// 	computeA(FlowVar::V);
	
// 	//Run the time steps.
// 	for(fV.timeId = 0; fV.timeId < 1; fV.timeId++)
// 	{
// 		computeConvectiveFlux(FlowVar::U);
// 		computeConvectiveFlux(FlowVar::V);

// 		computeDiffusiveFlux(FlowVar::U);
// 		computeDiffusiveFlux(FlowVar::V);

// 		computeb0(FlowVar::U);	//Compute the convection diffusion contribution to RHS of u
// 		computeb0(FlowVar::V);	//Compute the convection diffusion contribution to RHS of v

// 		//Run outer iterations
// 		while(m < 1)
// 		{
// 			solveXMom();
// 			solveYMom();

// 			fV.u = u_;
// 			fV.v = v_;
// 			m++;

// 		}



// 	}

// }

// void Driver::InitializeVelocity()
// {
// 	unsigned long Nx, Ny;
// 	Nx = G->Nx(), Ny = G->Ny();

// 	for(unsigned long i = 1; i < Nx-1; i++)
// 		for(unsigned long j = 1; j < Ny-1; j++){
// 			fV.u(i,j) = 0.0;//G->y(j) - 0.5;
// 			fV.v(i,j) = 0.0;//-(G->x(i) - 0.5);
// 		}
// }

// void Driver::InitializePressure()
// {
// 	unsigned long Nx, Ny;
// 	Nx = G->Nx(), Ny = G->Ny();
// 	for(unsigned long i = 1; i < Nx-1; i++)
// 		for(unsigned long j = 1; j < Ny-1; j++){
// 			fV.p(i,j) = 0.0;//pow(G->x(i),1) + pow(G->y(j),1);
// 		}
// }





// }

// void Driver::Test()
// {
// 	testConvectiveFlux();
// 	testDiffusiveFlux();
// 	testPressureGradient();
// }




// void initializeFlowVars(flowVars &fV, const Grid &G)
// {
// 	fV.timeId = 0;
// 	fV.time = 0.0;

// 	fV.u.setSize(G.Nx(),G.Ny());
// 	fV.v.setSize(G.Nx(),G.Ny());
// 	fV.p.setSize(G.Nx(),G.Ny());

// 	fV.uCn.setSize(G.Nx(),G.Ny()); fV.uCn_.setSize(G.Nx(),G.Ny());
// 	fV.uDn.setSize(G.Nx(),G.Ny());

// 	fV.vCn.setSize(G.Nx(),G.Ny()); fV.vCn_.setSize(G.Nx(),G.Ny());
// 	fV.vDn.setSize(G.Nx(),G.Ny());

// 	fV.gradPx.setSize(G.Nx(),G.Ny());
// 	fV.gradPy.setSize(G.Nx(),G.Ny());
// }