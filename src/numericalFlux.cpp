#include "cellState.h"
#include "grid.h"
#include <vector>
#include <cmath>
#include <numericalFlux.h>
#include "eulerFlux.h"

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
    vector<double> EFlux = EastFlux(grd,stencil);
    vector<double> EJAV = EastJamVisc(grd, stencil);
    vector<double> FluxAV(4, 0.0);

    FluxAV[0] = EFlux[0] - EJAV[0];
    FluxAV[1] = EFlux[1] - EJAV[1];
    FluxAV[2] = EFlux[2] - EJAV[2];
    FluxAV[3] = EFlux[3] - EJAV[3];
    return FluxAV;
}
vector<double> WestFlux_AV(grid &grd,  vector<cellState> &stencil){
    vector<double> WFlux = WestFlux(grd,stencil);
    vector<double> WJAV = WestJamVisc(grd, stencil);
    vector<double> FluxAV(4, 0.0);

    FluxAV[0] = WFlux[0] - WJAV[0];
    FluxAV[1] = WFlux[1] - WJAV[1];
    FluxAV[2] = WFlux[2] - WJAV[2];
    FluxAV[3] = WFlux[3] - WJAV[3];
    return FluxAV;
}
vector<double> NorthFlux_AV(grid &grd,  vector<cellState> &stencil){
    vector<double> NFlux = NorthFlux(grd, stencil);
    vector<double> NJAV = NorthJamVisc(grd, stencil);
    vector<double> FluxAV(4, 0.0);

    FluxAV[0] = NFlux[0] - NJAV[0];
    FluxAV[1] = NFlux[1] - NJAV[1];
    FluxAV[2] = NFlux[2] - NJAV[2];
    FluxAV[3] = NFlux[3] - NJAV[3];
    return FluxAV;
}
vector<double> SouthFlux_AV(grid &grd,  vector<cellState> &stencil){
    vector<double> SFlux = SouthFlux(grd, stencil);
    vector<double> SJAV = SouthJamVisc(grd, stencil);
    vector<double> FluxAV(4, 0.0);

    FluxAV[0] = SFlux[0] - SJAV[0];
    FluxAV[1] = SFlux[1] - SJAV[1];
    FluxAV[2] = SFlux[2] - SJAV[2];
    FluxAV[3] = SFlux[3] - SJAV[3];


    return FluxAV;
}

vector<double> AirfoilFlux(grid & grd, vector< vector<cellState> > &cellset, int i) {
	vector<double> FoilFlux(4, 0.0);
    vector<double> FSTAR(4,0.0);
    vector<double> GSTAR(4,0.0);
	cellState foil=cellset[i][0];
	double gamma=foil.gamma();
	double Pb = 3.0 / 2 * cellset[i][0].P() - 1.0 / 2 * cellset[i][1].P();// Linear extrapolation of pressure. See p 21 of stanford
    //Pb = cellset[i][0].P();//zeroth order approximation in space and time
    double rho=cellset[i][0].rho()*pow(cellset[i][0].P()/Pb, 1/gamma); //isentropic density change within cell

    double ubtang=foil.U()*grd.ySnorm[i][0] - foil.V()*grd.xSnorm[i][0];
    double ubx=ubtang*grd.ySnorm[i][0];
    double uby=ubtang*grd.xSnorm[i][0];

	FSTAR[0] = 0;
	FSTAR[1] = Pb*grd.xSnorm[i][0];
	FSTAR[2] = 0;
	FSTAR[3] = 0;

	GSTAR[0] = 0;
	GSTAR[1] = 0;
	GSTAR[2] = Pb*grd.ySnorm[i][0];
	GSTAR[3] = 0;

/* Removed to try the simple pressure normal flux
    FSTAR[0]= rho*ubx;
    FSTAR[1]= rho*pow(ubx,2)+Pb;
    FSTAR[2]= rho*ubx*uby;
    FSTAR[3]= ubx*(Pb*gamma/(gamma-1)+rho*(pow(ubx,2)+ pow(uby, 2))/2);

    GSTAR[0]= rho*uby;
    GSTAR[1]= rho*ubx*uby;
    GSTAR[2]= rho*pow(uby,2)+Pb;
	GSTAR[3]= uby*(Pb*gamma / (gamma - 1) + rho*(pow(ubx, 2) + pow(uby, 2)) / 2);

	*/

    FoilFlux[0]=FSTAR[0]*grd.xSdeltas[i][0]+GSTAR[0]*grd.ySdeltas[i][0];
    FoilFlux[1]=FSTAR[1]*grd.xSdeltas[i][0]+GSTAR[1]*grd.ySdeltas[i][0];
    FoilFlux[2]=FSTAR[2]*grd.xSdeltas[i][0]+GSTAR[2]*grd.ySdeltas[i][0];
    FoilFlux[3]=FSTAR[3]*grd.xSdeltas[i][0]+GSTAR[3]*grd.ySdeltas[i][0];

//    FoilFlux[0] = foil.F1()*grd.xSdeltas[i][0] + foil.G1()*grd.ySdeltas[i][0];
//    FoilFlux[1] = foil.F2()*grd.xSdeltas[i][0] + foil.G2()*grd.ySdeltas[i][0];
//    FoilFlux[2] = foil.F3()*grd.xSdeltas[i][0] + foil.G3()*grd.ySdeltas[i][0];
//    FoilFlux[3] = foil.F4()*grd.xSdeltas[i][0] + foil.G4()*grd.ySdeltas[i][0]; // THese should definitely be delta S's

//	FoilFlux[0] = 0;
//	FoilFlux[1] = Pboundary*grd.xSnorm[i][0];
//	FoilFlux[2] = Pboundary*grd.ySnorm[i][0];
//	FoilFlux[3] = 0;

//    vector<double> test(4,0.0);
//    test[0]=ubx;
//    test[1]=uby;
//    test[2]=ubtang;
//    test[3]=Pb;
	return FoilFlux;
}

vector<double> InletOutletFlux(grid &grd,  vector< vector<cellState> > &cellset, int i){
    //Reference values (we should maybe find a less janky way of getting these in here)
    double M_ref=0.8;
    double rho_ref=1.293; // Kg/m^3
    double P_ref=101325;// Pa
    double c_ref=sqrt(1.4*P_ref/rho_ref);
    double R=8314./28.966; // air gas constant J/Kg-K
    double cp=1005.; // J/Kg-K
    double cv=cp-R;
    double gamma=cp/cv;
    double rhoE_ref=P_ref/(gamma-1)+0.5*rho_ref*pow(M_ref*c_ref,2); //internal energy KJ/m^3,

    double Rplus, Rminus, ubnorm,ubtang, ubx, uby, ceebee, rho, P;
    vector<double> FSTAR(4,0.0);
    vector<double> GSTAR(4,0.0);
    vector<double> netFlux(4,0.0);


    int ctrlast=grd.M-2; // Last is actually M-2, unless we give cellset a ghost cell
	int normlast = grd.M - 1; // The normal on the boundary is what matters here

    double normal_speed=cellset[i][ctrlast].U()*grd.xSnorm[i][normlast] + cellset[i][ctrlast].V()*grd.ySnorm[i][normlast];

    if(normal_speed < 0){ //inlet
//        Rplus= M_ref*c_ref*grd.xSnorm[i][normlast]	 +2*cellset[i][ctrlast].C()/(gamma-1);
//        Rminus= cellset[i][ctrlast].U()*grd.xSnorm[i][normlast]+ cellset[i][ctrlast].V()*grd.ySnorm[i][normlast]-2*cellset[i][ctrlast].C()/(gamma-1); //physical inlet velcity dotted with the normal vector
//        ubnorm=(Rplus+Rminus)/2; //normal component of velocity
//        ubtang=cellset[i][ctrlast].U()*grd.ySnorm[i][normlast]-cellset[i][ctrlast].V()*grd.xSnorm[i][normlast]; //tangential component. Geometry looks good.
//        ubx = (ubnorm-ubtang)*grd.xSnorm[i][normlast]; //check geometrically
//        uby = (ubnorm+ubtang)*grd.ySnorm[i][normlast]; //check geometrically
//        ceebee = (Rplus - Rminus)*(gamma - 1) / 4;
		P=P_ref;
		ubx=M_ref*c_ref;
		uby=0;

//        cellState freestream=cellState(rho_ref, rho_ref*M_ref*c_ref, 0, rhoE_ref, gamma, cv, i,-1);
//        vector<cellState> stencil(6);
//        stencil[3]=freestream;
//        stencil[2]=cellset[i][normlast];
//        return GenericFlux(grd, stencil);
    }
    else{//outlet
        Rminus=normal_speed+2*cellset[i][ctrlast].C()/(gamma-1);
        Rplus= M_ref*c_ref*grd.xSnorm[i][normlast]+2*c_ref/(gamma-1); //physical outlet velcity dotted with the normal vector (out outlet velocity is only in the x-driection)
        ubnorm=(Rplus+Rminus)/2; //normal component of velocity
        ubtang=cellset[i][ctrlast].U()*grd.ySnorm[i][normlast]-cellset[i][ctrlast].V()*grd.xSnorm[i][normlast]; //tangential component. Geometry looks good.
        ubx = (ubnorm-ubtang)*grd.xSnorm[i][normlast];
        uby = (ubnorm+ubtang)*grd.ySnorm[i][normlast];
        ubx=M_ref*c_ref; //X-speed is at the mach speed, ignoring ferrante's wild notions of consistant tangential velocity
        uby=0;
        ceebee = (Rminus - Rplus)*(gamma - 1) / 4;
    }

    //rho = pow(pow(ceebee, 2)*pow(rho_ref, gamma) / (gamma*P_ref), gamma - 1);
    //P = pow(ceebee, 2)*rho / gamma;
    P=P_ref;
    rho=rho_ref;

    FSTAR[0]= rho*ubx;
    FSTAR[1]= rho*pow(ubx,2)+P;
    FSTAR[2]= rho*ubx*uby;
    FSTAR[3]= ubx*(P*gamma/(gamma-1)+rho*(pow(ubx,2)+ pow(uby, 2))/2);

    GSTAR[0]= rho*uby;
    GSTAR[1]= rho*ubx*uby;
    GSTAR[2]= rho*pow(uby,2)+P;
	GSTAR[3]= uby*(P*gamma / (gamma - 1) + rho*(pow(ubx, 2) + pow(uby, 2)) / 2);

    netFlux[0]=FSTAR[0]*grd.xSdeltas[i][normlast]+GSTAR[0]*grd.ySdeltas[i][normlast];
    netFlux[1]=FSTAR[1]*grd.xSdeltas[i][normlast]+GSTAR[1]*grd.ySdeltas[i][normlast];
    netFlux[2]=FSTAR[2]*grd.xSdeltas[i][normlast]+GSTAR[2]*grd.ySdeltas[i][normlast];
    netFlux[3]=FSTAR[3]*grd.xSdeltas[i][normlast]+GSTAR[3]*grd.ySdeltas[i][normlast];

//    vector<double> test(4,0.0);
//    test[0]=FSTAR[0]*grd.xSnorm[i][normlast]+GSTAR[0]*grd.ySnorm[i][normlast];
//    test[1]=FSTAR[1]*grd.xSnorm[i][normlast]+GSTAR[1]*grd.ySnorm[i][normlast];
//    test[2]=FSTAR[2]*grd.xSnorm[i][normlast]+GSTAR[2]*grd.ySnorm[i][normlast];
//    test[3]=FSTAR[3]*grd.xSnorm[i][normlast]+GSTAR[3]*grd.ySnorm[i][normlast];
//    test[0]=(normal_speed/abs(normal_speed))*sqrt(ubx*ubx+uby*uby);
//    test[1]=ceebee;
//    test[2]=P;
//    test[3]=(P/rho*gamma / (gamma - 1) + (pow(ubx, 2) + pow(uby, 2)) / 2);

    return netFlux;
}
