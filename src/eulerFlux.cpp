#include "cellState.h"
#include "grid.h"
#include <vector>
#include <iostream>
#include "eulerFlux.h"

using namespace std;
vector<double> GenericFlux(grid &grd,  vector<cellState> &stencil){
    double FSTAR1, FSTAR2, FSTAR3, FSTAR4;
    double GSTAR1, GSTAR2, GSTAR3, GSTAR4;

    FSTAR1 = 0.5*(stencil[2].F1() + stencil[3].F1());
	FSTAR2 = 0.5*(stencil[2].F2() + stencil[3].F2());
	FSTAR3 = 0.5*(stencil[2].F3() + stencil[3].F3());
	FSTAR4 = 0.5*(stencil[2].F4() + stencil[3].F4());
	GSTAR1 = 0.5*(stencil[2].G1() + stencil[3].G1());
	GSTAR2 = 0.5*(stencil[2].G2() + stencil[3].G2());
	GSTAR3 = 0.5*(stencil[2].G3() + stencil[3].G3());
	GSTAR4 = 0.5*(stencil[2].G4() + stencil[3].G4());

	int I=stencil[3].i();
    int J=stencil[3].j();

	std::vector<double> GENERIC_FLUX(4, 0.0);
	//Perpendicular flux to the boundary normal is <F, G> dot <nx, ny>, but eventually, deltaS must be added too
	GENERIC_FLUX[0] = FSTAR1*grd.xWdeltas[I][J] + GSTAR1*grd.yWdeltas[I][J];
	GENERIC_FLUX[1] = FSTAR2*grd.xWdeltas[I][J] + GSTAR2*grd.yWdeltas[I][J];
	GENERIC_FLUX[2] = FSTAR3*grd.xWdeltas[I][J] + GSTAR3*grd.yWdeltas[I][J];
	GENERIC_FLUX[3] = FSTAR4*grd.xWdeltas[I][J] + GSTAR4*grd.yWdeltas[I][J]; // THese should definitely be delta S's
	return GENERIC_FLUX;

}

vector<double> EastFlux(grid &grd, vector<cellState> &stencil){
    vector<cellState> substencil(6);
    for(int i=0; i<6; i++){
        substencil[i]=stencil[i+1];
    }
    return GenericFlux(grd, substencil);
}
vector<double> WestFlux(grid &grd, vector<cellState> &stencil){
    vector<cellState> substencil(6);
    for(int i=0; i<6; i++){
        substencil[i]=stencil[i];
    }
    return GenericFlux(grd, substencil);
}
vector<double> NorthFlux(grid &grd, vector<cellState> &stencil){
    vector<cellState> substencil(6);
    for(int i=0; i<6; i++){
        substencil[i]=stencil[i+1];
    }
    return GenericFlux(grd, substencil);
}
vector<double> SouthFlux(grid &grd, vector<cellState> &stencil){
    vector<cellState> substencil(6);
    for(int i=0; i<6; i++){
        substencil[i]=stencil[i];
    }
    return GenericFlux(grd, substencil);
}

