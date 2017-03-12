#include "cellState.h"
#include "grid.h"
#include <vector>
#include <cmath>
#include <numericalFlux.h>

using namespace std;

double nu_max(vector<cellState> stencil){
    double nu_max = 0.0;
	double nu_test;

	//This loop to find the maximum nu iterates through k seeing if nu[i+k] is greater than any previous contenders of nu
	for (int k = 0; k<=3; k++) {
		//the expression for nu with (i+k) in place of i to allow the loop to examine nu(i-1), nu(i), nu(i+1), and nu(i+2)
		nu_test = abs(stencil[2+k].P() - 2*stencil[1+k].P()  + stencil[k].P());
		nu_test /= stencil[2+k].P() + 2*stencil[1+k].P()  + stencil[k].P();
		if (nu_test > nu_max){
			nu_max = nu_test;
		}
	}
	return nu_max;
}

vector<double> GenericJamesonViscocity(grid &grd, vector<cellState> &stencil) {
    //Generic jameson viscocity focused upwind on the cell in stencil[2] in which stencil is a vector of 6 cells

    double alpha2=1.0/4;
    double alpha4=1.0/256;

    //Calculate the epsilon values at i+1/2,j
    double epsilon2, epsilon4;
    // prepare the (i+1/2,j) values needed in the epsilon terms
    double u_avg = (stencil[2].U() + stencil[3].U()) / 2;
    double v_avg = (stencil[2].V() + stencil[3].V()) / 2;
    double c_avg = (stencil[2].C() + stencil[3].C()) / 2;

    // The dot product term in the epsilon terms
    double delta_sx=grd.xWnorm[stencil[3].i()][stencil[3].j()]*grd.yWside[stencil[3].i()][stencil[3].j()];
    double delta_sy=grd.yWnorm[stencil[3].i()][stencil[3].j()]*grd.xWside[stencil[3].i()][stencil[3].j()];
    double u_vec_dot_s_vec = u_avg*delta_sx + v_avg*delta_sy;

    // Solve for epsilon 2 and 4
    epsilon2 = 0.5*alpha2*(u_vec_dot_s_vec + c_avg)*nu_max(stencil);
    epsilon4 = 0.5*alpha4*(u_vec_dot_s_vec + c_avg) - epsilon2;
    if (epsilon4 < 0.0) {
        epsilon4 = 0.0;
    }
    vector<double> D_AV(4, 0.0), dUds(4, 0.0), d3Uds3(4, 0.0);

    dUds[0]=stencil[3].rho() - stencil[2].rho();
    dUds[1]=stencil[3].rhoU() - stencil[2].rhoU();
    dUds[2]=stencil[3].rhoV() - stencil[2].rhoV();
    dUds[3]=stencil[3].rhoE() - stencil[2].rhoE();

    d3Uds3[0]=stencil[4].rho() - 3*stencil[3].rho() + 3*stencil[2].rho() - stencil[1].rho();
    d3Uds3[1]=stencil[4].rhoU() - 3*stencil[3].rhoU() + 3*stencil[2].rhoU() - stencil[1].rhoU();
    d3Uds3[2]=stencil[4].rhoV() - 3*stencil[3].rhoV() + 3*stencil[2].rhoV() - stencil[1].rhoV();
    d3Uds3[3]=stencil[4].rhoE() - 3*stencil[3].rhoE() + 3*stencil[2].rhoE() - stencil[1].rhoE();

    // The Jamison Artificial Viscosity
    D_AV[0] = epsilon2*dUds[0] - epsilon4*d3Uds3[0];
    D_AV[1] = epsilon2*dUds[1] - epsilon4*d3Uds3[1];
    D_AV[2] = epsilon2*dUds[2] - epsilon4*d3Uds3[2];
    D_AV[3] = epsilon2*dUds[3] - epsilon4*d3Uds3[3];
    return D_AV;
}

vector<double> EastJamVisc(grid &grd,  vector<cellState> &stencil) {
    vector<cellState> substencil(6);
    for(int i=0; i<6; i++){
        substencil[i]=stencil[i+1];
    }
    return GenericJamesonViscocity(grd, substencil);
    //return GenericJamVisc(grd, cell i-1, cell i, cell i+1, cell i+2);
}
vector<double> WestJamVisc(grid &grd,  vector<cellState> &stencil) {
    vector<cellState> substencil(6);
    for(int i=0; i<6; i++){
        substencil[i]=stencil[i];
    }
    return GenericJamesonViscocity(grd, substencil);
}
vector<double> NorthJamVisc(grid &grd,  vector<cellState> &stencil) {
    vector<cellState> substencil(6);
    for(int i=0; i<6; i++){
        substencil[i]=stencil[i+1];
    }
    return GenericJamesonViscocity(grd, substencil);
}
vector<double> SouthJamVisc(grid &grd,  vector<cellState> &stencil) {
    vector<cellState> substencil(6);
    for(int i=0; i<6; i++){
        substencil[i]=stencil[i];
    }
    return GenericJamesonViscocity(grd, substencil);
    // is equivilent to return GenericJamesonViscocity(grd, cellset, i, j-1, 0, 1)
}

vector<double> EastFlux_AV(grid &grd,  vector<cellState> &stencil){

    vector<double> EFlux = EastFlux(grd,cellset,i,j);
    vector<double> EJAV = EastJamVisc(grd, cellset,i,j);
    vector<double> FluxAV(4, 0.0);

    FluxAV[0] = EFlux[0] - EJAV[0];
    FluxAV[1] = EFlux[1] - EJAV[1];
    FluxAV[2] = EFlux[2] - EJAV[2];
    FluxAV[3] = EFlux[3] - EJAV[3];


    return FluxAV;
}
vector<double> WestFlux_AV(grid &grd, vector< vector<cellState> > &cellset,int i, int j){

    vector<double> WFlux = WestFlux(grd,cellset,i,j);
    vector<double> WJAV = WestJamVisc(grd, cellset,i,j);
    vector<double> FluxAV(4, 0.0);

    FluxAV[0] = WFlux[0] - WJAV[0];
    FluxAV[1] = WFlux[1] - WJAV[1];
    FluxAV[2] = WFlux[2] - WJAV[2];
    FluxAV[3] = WFlux[3] - WJAV[3];


    return FluxAV;
}
vector<double> NorthFlux_AV(grid &grd, vector< vector<cellState> > &cellset,int i, int j){

    vector<double> NFlux = NorthFlux(grd,cellset,i,j);
    vector<double> NJAV = NorthJamVisc(grd, cellset,i,j);
    vector<double> FluxAV(4, 0.0);

    FluxAV[0] = NFlux[0] - NJAV[0];
    FluxAV[1] = NFlux[1] - NJAV[1];
    FluxAV[2] = NFlux[2] - NJAV[2];
    FluxAV[3] = NFlux[3] - NJAV[3];


    return FluxAV;
}
vector<double> SouthFlux_AV(grid &grd, vector< vector<cellState> > &cellset,int i, int j){

    vector<double> SFlux = SouthFlux(grd,cellset,i,j);
    vector<double> SJAV = SouthJamVisc(grd, cellset,i,j);
    vector<double> FluxAV(4, 0.0);

    FluxAV[0] = SFlux[0] - SJAV[0];
    FluxAV[1] = SFlux[1] - SJAV[1];
    FluxAV[2] = SFlux[2] - SJAV[2];
    FluxAV[3] = SFlux[3] - SJAV[3];


    return FluxAV;
}

vector<double> Residuals(grid &grd, vector< vector<cellState> > &cellset,int i, int j){

    std::vector<double> RESIDUALS(4, 0.0);
    std::vector<double> NFAV(4, 0.0);
    std::vector<double> SFAV(4, 0.0);
    std::vector<double> EFAV(4, 0.0);
    std::vector<double> WFAV(4, 0.0);

    NFAV = NorthFlux_AV(grd,cellset,i,j);
    SFAV = SouthFlux(grd,cellset,i,j);
    EFAV = EastFlux_AV(grd,cellset,i,j) ;
    WFAV = WestFlux_AV(grd,cellset,i,j);

    RESIDUALS[0] = NFAV[0] - SFAV[0] + EFAV[0] - WFAV[0];
    RESIDUALS[1] = NFAV[1] - SFAV[1] + EFAV[1] - WFAV[1];
    RESIDUALS[2] = NFAV[2] - SFAV[2] + EFAV[2] - WFAV[2];
    RESIDUALS[3] = NFAV[3] - SFAV[3] + EFAV[3] - WFAV[3];

    return RESIDUALS;

}

double Tau(grid &grd, vector< vector<cellState> > &cellset, int i, int j, double CFL){
    double denominatorI = (cellset[i][j].U()+cellset[i][j].C())*grd.xInorm[i][j] + (cellset[i][j].V()+cellset[i][j].C())*grd.yInorm[i][j];
    double denominatorJ = (cellset[i][j].U()+cellset[i][j].C())*grd.xJnorm[i][j] + (cellset[i][j].V()+cellset[i][j].C())*grd.yJnorm[i][j];
    return CFL*grd.area[i][j]/(abs(denominatorI)+abs(denominatorJ)); // ------------- Consider making Tau 0.75 of this if it matters
}

vector<double> AlphaRK(){

    vector<double> Alpha;

    Alpha.push_back(0.125);
    Alpha.push_back(0.306);
    Alpha.push_back(0.587);
    Alpha.push_back(1.000);

    return Alpha;
}

//Define AlphaRK, Define Tau, Is there a way to do what I am trying to do with TempCell = CellSet?

vector< vector<cellState> > RK4(grid &grd, vector< vector<cellState> > &cellset, double CFL){

    std::vector<double> U_temp(4, 0.0);
    std::vector<double> RESIDUALS(4, 0.0);
    vector< vector<cellState> > cellsetPlus = cellset;
	vector< vector<cellState> > cellsetPrev(grd.N - 1, std::vector<cellState>(grd.M - 1));
    vector<double> alphaRK = AlphaRK();
    //loop through the whole thing 4 times! Nk*Ni*Nj, this way each pseudo timestep residual (R1,R2, etc) is based on fluxes from neighbors on that pseudotime. Otherwise there is the inclusion of fluxes that are old since Fstar is 1/2(fi + fi+1).
    for (int k = 0; k<4; k++) {
        cellsetPrev = cellsetPlus;
        for(int i=0+5; i<grd.N-1-5; i++) {
            for(int j=0+5; j<grd.M-1-5; j++) {
                RESIDUALS = Residuals(grd, cellsetPrev, i, j);

                U_temp[0] = cellset[i][j].rho() - Tau(grd,cellset,i,j,CFL) * alphaRK[k] * RESIDUALS[0];
                U_temp[1] = cellset[i][j].rhoU() - Tau(grd,cellset,i,j,CFL) * alphaRK[k] * RESIDUALS[1];
                U_temp[2] = cellset[i][j].rhoV() - Tau(grd,cellset,i,j,CFL) * alphaRK[k] * RESIDUALS[2];
                U_temp[3] = cellset[i][j].rhoE() - Tau(grd,cellset,i,j,CFL) * alphaRK[k] * RESIDUALS[3];

                cellsetPlus[i][j] = cellState(U_temp[0], U_temp[1], U_temp[2], U_temp[3], cellset[i][j].gamma(), cellset[i][j].cv(), i, j);
            }
        }
    }
    return cellsetPlus;
}
