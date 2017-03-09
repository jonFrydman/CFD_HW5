// Mike Wennerstrom - CFD 543 - Numerical Flux Function

// This function gives back the numerical flux. Jon, things will be redefined over and over if this loops. Is there a way around this? Maybe it doesn't matter.

// INPUTS:
// Cellset: a 2D vector of the cell states
// GeoSet: a 2D vector of geometric grid data
// i and j: integers to allow looping

// OUTPUTS: a 1D vector with each numerical flux. F[0]=rhoU, ... , F[3] = rho U H

#include "referenceParameters.h"
#include "cellState.h"
#include "grid.h"
#include<iostream>
#include <vector>
using namespace std;

vector<double> NumericalFlux(int i, int j) {

	double FSTAR1, FSTAR2, FSTAR3, FSTAR4;
    double GSTAR1, GSTAR2, GSTAR3, GSTAR4;
	double alpha2 = 0.25;
	double alpha4 = 1.0 / 256;

    //EastSide Fluxes on cell i,j
	FSTAR1 = 0.5*(cellset[i][j].F1() + cellset[i + 1][j].F1());
	FSTAR2 = 0.5*(cellset[i][j].F2() + cellset[i + 1][j].F2());
	FSTAR3 = 0.5*(cellset[i][j].F3() + cellset[i + 1][j].F3());
	FSTAR4 = 0.5*(cellset[i][j].F4() + cellset[i + 1][j].F4());
	GSTAR1 = 0.5*(cellset[i][j].G1() + cellset[i + 1][j].G1());
	GSTAR2 = 0.5*(cellset[i][j].G2() + cellset[i + 1][j].G2());
	GSTAR3 = 0.5*(cellset[i][j].G3() + cellset[i + 1][j].G3());
	GSTAR4 = 0.5*(cellset[i][j].G4() + cellset[i + 1][j].G4());

	//find max of (nu(i-1), nu(i), nu(i+1), nu(i+2))
	double nu_max = 0.0;
	double nu_test;

	//This loop to find the maximum nu iterates through k seeing if nu[i+k] is greater than any previous contenders of nu
	for (int k = -1; k<2; k++) {
		//the expression for nu with (i+k) in place of i to allow the loop to examine nu(i-1), nu(i), nu(i+1), and nu(i+2)
		nu_test = abs((cellset[(i+k)+ 1][j].P() - 2 * cellset[i+k][j].P() + cellset[(i+k) - 1][j].P()) / (cellset[(i+k) + 1][j].P() + 2 * cellset[i+k][j].P() + cellset[(i+k) - 1][j].P()));
		if (nu_test > nu_max) {
			nu_max = nu_test;
		}
	}
    //Calculate the epsilon values at i+1/2,j
	double epsilon2, epsilon4;
	// prepare the (i+1/2,j) values needed in the epsilon terms
	double u_avg = (cellset[i][j].U() + cellset[i + 1][j].U()) / 2;
	double v_avg = (cellset[i][j].V() + cellset[i + 1][j].V()) / 2;
	double c_avg = (cellset[i][j].C() + cellset[i + 1][j].C()) / 2;

	// The dot product term in the epsilon terms
	double u_vec_dot_s_vec = u_avg*grd.xWnorm[i + 1][j] + v_avg*grd.yWnorm[i + 1][j];

	// Solve for epsilon 2 and 4
	epsilon2 = 0.5*alpha2*(u_vec_dot_s_vec + c_avg)*nu_max;
	epsilon4 = 0.5*alpha4*(u_vec_dot_s_vec + c_avg) - epsilon2;
	if (epsilon4 < 0.0) {
		epsilon4 = 0.0;
	}

    vector<double> NUM_FLUX(4, 0.0);
	// The Jamison Artificial Viscosity
	double D_AV_1 = epsilon2*(cellset[i + 1][j].rho() - cellset[i][j].rho()) - epsilon4*(cellset[i + 2][j].rho() - 3 * cellset[i + 1][j].rho() + 3 * cellset[i][j].rho() - cellset[i - 1][j].rho());
	double D_AV_2 = epsilon2*(cellset[i + 1][j].rhoU() - cellset[i][j].rhoU()) - epsilon4*(cellset[i + 2][j].rhoU() - 3 * cellset[i + 1][j].rhoU() + 3 * cellset[i][j].rhoU() - cellset[i - 1][j].rhoU());
	double D_AV_3 = epsilon2*(cellset[i + 1][j].rhoV() - cellset[i][j].rhoV()) - epsilon4*(cellset[i + 2][j].rhoV() - 3 * cellset[i + 1][j].rhoV() + 3 * cellset[i][j].rhoV() - cellset[i - 1][j].rhoV());
	double D_AV_4 = epsilon2*(cellset[i + 1][j].rhoE() - cellset[i][j].rhoE()) - epsilon4*(cellset[i + 2][j].rhoE() - 3 * cellset[i + 1][j].rhoE() + 3 * cellset[i][j].rhoE() - cellset[i - 1][j].rhoE());

	// EastSide Flux on cell i,j
	NUM_FLUX[0] = FSTAR1*grd.xWnorm[i + 1][j] + GSTAR1*grd.yWnorm[i+1][j];
	NUM_FLUX[1] = FSTAR2*grd.xWnorm[i + 1][j] + GSTAR2*grd.yWnorm[i + 1][j];
	NUM_FLUX[2] = FSTAR3*grd.xWnorm[i + 1][j] + GSTAR3*grd.yWnorm[i + 1][j];
	NUM_FLUX[3] = FSTAR4*grd.xWnorm[i + 1][j] + GSTAR4*grd.yWnorm[i + 1][j];

	return NUM_FLUX;
}
