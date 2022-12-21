#include "Base.h"

#ifndef CONVECTION_DIFFUSION_H
#define CONVECTION_DIFFUSION_H

void computeConvectiveXFlux(Vector &Fc, const Vector &u, const Vector &v,  const Grid &G);
void computeConvectiveYFlux(Vector &Fc, const Vector &u, const Vector &v,  const Grid &G);
void computeDiffusiveFlux(Vector &Fd, const Vector &phi, const Grid &G);
void computeGradient(Vector &gradphix, Vector&gradphiy, const Vector &phi, const Grid &G);

void testFlux(const Grid &G);






#endif