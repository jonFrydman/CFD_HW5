#include "cellState.h"
#include "grid.h"
#include <vector>
#include <cmath>
#include <numericalFlux.h>

using namespace std;

vector<double> GenericFlux(grid &grd, vector< vector<cellState> > &cellset, int i, int j, int delta_i, int delta_j){
    double FSTAR1, FSTAR2, FSTAR3, FSTAR4;
    double GSTAR1, GSTAR2, GSTAR3, GSTAR4;

    FSTAR1 = 0.5*(cellset[i][j].F1() + cellset[i + delta_i][j+delta_j].F1());
	FSTAR2 = 0.5*(cellset[i][j].F2() + cellset[i + delta_i][j+delta_j].F2());
	FSTAR3 = 0.5*(cellset[i][j].F3() + cellset[i + delta_i][j+delta_j].F3());
	FSTAR4 = 0.5*(cellset[i][j].F4() + cellset[i + delta_i][j+delta_j].F4());
	GSTAR1 = 0.5*(cellset[i][j].G1() + cellset[i + delta_i][j+delta_j].G1());
	GSTAR2 = 0.5*(cellset[i][j].G2() + cellset[i + delta_i][j+delta_j].G2());
	GSTAR3 = 0.5*(cellset[i][j].G3() + cellset[i + delta_i][j+delta_j].G3());
	GSTAR4 = 0.5*(cellset[i][j].G4() + cellset[i + delta_i][j+delta_j].G4());

	std::vector<double> GENERIC_FLUX(4, 0.0);

	GENERIC_FLUX[0] = FSTAR1*grd.xWnorm[i+delta_i][j+delta_j] + GSTAR1*grd.yWnorm[i+delta_i][j+delta_j];
	GENERIC_FLUX[1] = FSTAR2*grd.xWnorm[i+delta_i][j+delta_j] + GSTAR2*grd.yWnorm[i+delta_i][j+delta_j];
	GENERIC_FLUX[2] = FSTAR3*grd.xWnorm[i+delta_i][j+delta_j] + GSTAR3*grd.yWnorm[i+delta_i][j+delta_j];
	GENERIC_FLUX[3] = FSTAR4*grd.xWnorm[i+delta_i][j+delta_j] + GSTAR4*grd.yWnorm[i+delta_i][j+delta_j];

	return GENERIC_FLUX;
}

vector<double> EastFlux(grid &grd, vector< vector<cellState> > &cellset, int i, int j){
    return GenericFlux(grd, cellset, i, j, 1, 0);
}
vector<double> WestFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j){
    return GenericFlux(grd, cellset, i, j, -1, 0);
}
vector<double> NorthFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j){
    return GenericFlux(grd, cellset, i, j, 0, 1);
}
vector<double> SouthFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j){
    return GenericFlux(grd, cellset, i, j, 0, -1);
}

vector<double> JamesonViscocity(grid &grd, vector< vector<cellState> > &cellset,int i, int j) {
	//find max of (nu(i-1), nu(i), nu(i+1), nu(i+2))
	double nu_max = 0.0;
	double nu_test;

	//This loop to find the maximum nu iterates through k seeing if nu[i+k] is greater than any previous contenders of nu
	for (int k = -1; k<2; k++) {
		//the expression for nu with (i+k) in place of i to allow the loop to examine nu(i-1), nu(i), nu(i+1), and nu(i+2)
		nu_test = abs((cellset[(i+k)+ 1][j].P() - 2 * cellset[i+k][j].P() + cellset[(i+k) - 1][j].P()) / (cellset[(i+k) + 1][j].P() + 2 * cellset[i+k][j].P() + cellset[(i+k) - 1][j].P()));
		if (nu_test > nu_max){
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
    vector<double> D_AV(4, 0.0);
	// The Jamison Artificial Viscosity
	D_AV[0] = epsilon2*(cellset[i + 1][j].rho() - cellset[i][j].rho()) - epsilon4*(cellset[i + 2][j].rho() - 3 * cellset[i + 1][j].rho() + 3 * cellset[i][j].rho() - cellset[i - 1][j].rho());
	D_AV[1] = epsilon2*(cellset[i + 1][j].rhoU() - cellset[i][j].rhoU()) - epsilon4*(cellset[i + 2][j].rhoU() - 3 * cellset[i + 1][j].rhoU() + 3 * cellset[i][j].rhoU() - cellset[i - 1][j].rhoU());
	D_AV[2] = epsilon2*(cellset[i + 1][j].rhoV() - cellset[i][j].rhoV()) - epsilon4*(cellset[i + 2][j].rhoV() - 3 * cellset[i + 1][j].rhoV() + 3 * cellset[i][j].rhoV() - cellset[i - 1][j].rhoV());
	D_AV[3] = epsilon2*(cellset[i + 1][j].rhoE() - cellset[i][j].rhoE()) - epsilon4*(cellset[i + 2][j].rhoE() - 3 * cellset[i + 1][j].rhoE() + 3 * cellset[i][j].rhoE() - cellset[i - 1][j].rhoE());
	return D_AV;
}

//void IFT_IICD(int index, double dX, double dT, std::vector<double> & uT, std::vector<double> & uTplus);
//void IBT_2CD_I(int Xelements, std::vector<double>& CoefficientMinus, std::vector<double>& CoefficientPlus, std::vector<double>& Coefficient, std::vector<double>& Uplus, std::vector<double>& U);

