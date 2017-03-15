#include "integration.h"
#include "cellState.h"
#include "grid.h"
#include <vector>
#include <cmath>
#include <numericalFlux.h>

vector<double> Residuals(grid &grd,  vector< vector<cellState> > &cellset, int i, int j){

    std::vector<double> RESIDUALS(4, 0.0);
    std::vector<double> NFAV(4, 0.0);
    std::vector<double> SFAV(4, 0.0);
    std::vector<double> EFAV(4, 0.0);
    std::vector<double> WFAV(4, 0.0);

	std::cout << i << ' ' << j <<endl;

    vector<cellState> ns_stencil=stencilNS(grd, cellset, i, j);
    vector<cellState> ew_stencil=stencilEW(grd, cellset, i, j);
    if(j>=2 && j<grd.M-4){ //if away from the airfoil and the airfoil boundary, use artificial viscocity
        NFAV = NorthFlux_AV(grd, ns_stencil);
        SFAV = SouthFlux_AV(grd, ns_stencil);
    }
    else if(j==0){ //at airfoil boundary airfoil case, no AV
        NFAV = NorthFlux(grd, ns_stencil);
        SFAV = AirfoilFlux(grd, cellset, i);
    }
    else if(j==grd.M-2){

        NFAV = InletOutletFlux(grd, cellset, i);
        SFAV = SouthFlux(grd, ns_stencil);
    
	}
    else{ //near boundaries assume AV=0
        NFAV = NorthFlux(grd, ns_stencil);
        SFAV = SouthFlux(grd, ns_stencil);
    }
    //AV always on in the East/West direction
    EFAV = EastFlux_AV(grd, ew_stencil) ;
    WFAV = WestFlux_AV(grd, ew_stencil);

    RESIDUALS[0] = NFAV[0] - SFAV[0] + EFAV[0] - WFAV[0];
    RESIDUALS[1] = NFAV[1] - SFAV[1] + EFAV[1] - WFAV[1];
    RESIDUALS[2] = NFAV[2] - SFAV[2] + EFAV[2] - WFAV[2];
    RESIDUALS[3] = NFAV[3] - SFAV[3] + EFAV[3] - WFAV[3];
    return RESIDUALS;
}

double Tau(grid &grd, cellState cell, double CFL){
    int i= cell.i();
    int j= cell.j();
    double denominatorI = (cell.U()+cell.C())*grd.xInorm[i][j] + (cell.V()+cell.C())*grd.yInorm[i][j]; //Yikes, Inorm isn't calculated for 127??
    double denominatorJ = (cell.U()+cell.C())*grd.xJnorm[i][j] + (cell.V()+cell.C())*grd.yJnorm[i][j];
    
	double TAU = CFL*grd.area[i][j] / (abs(denominatorI) + abs(denominatorJ));
	
	return TAU;
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
	double TauVal;

    //loop through the whole thing 4 times! Nk*Ni*Nj, this way each pseudo timestep residual (R1,R2, etc) is based on fluxes from neighbors on that pseudotime. Otherwise there is the inclusion of fluxes that are old since Fstar is 1/2(fi + fi+1).
    for (int k = 0; k<4; k++) {
        cellsetPrev = cellsetPlus;
            for(int j=0; j<grd.M-1; j++) { // not sure why i was looped first. 
				for (int i = 0; i<grd.N - 1; i++) {
				RESIDUALS = Residuals(grd, cellsetPrev, i, j);

				TauVal = Tau(grd, cellset[i][j], CFL);

                U_temp[0] = cellset[i][j].rho() - TauVal * alphaRK[k] * RESIDUALS[0];
                U_temp[1] = cellset[i][j].rhoU() - TauVal * alphaRK[k] * RESIDUALS[1];
                U_temp[2] = cellset[i][j].rhoV() - TauVal * alphaRK[k] * RESIDUALS[2];
                U_temp[3] = cellset[i][j].rhoE() - TauVal * alphaRK[k] * RESIDUALS[3];

                cellsetPlus[i][j] = cellState(U_temp[0], U_temp[1], U_temp[2], U_temp[3], cellset[i][j].gamma(), cellset[i][j].cv(), i, j);
            }
        }
    }
    return cellsetPlus;
}

vector<cellState> stencilEW(grid & grd, vector< vector<cellState> > & cellset, int i, int j){
    vector<cellState> stencil(7);

    if(i+3<cellset.size() && i-3 >=0){//Interior stencil
        for(int n=0; n<7; n++){
            stencil[n]=cellset[i+n-3][j];
        }
    }
    else if(i+3 >= cellset.size()){ //Stencil focused near the I=IMAX border
        int cnt=i-3;
        int n=0;
        while(n<7){
            stencil[n++]=cellset[cnt++][j]; //I'm not sure n++ accounts for stencil[0] at i=127 since stencil[4] = cellset[128][j] => nan Since cell 128 is a ghost cell
            if(cnt>=cellset.size()-1){ // ************** Changed this from cnt>=cellset.size() so that when cnt=127, and cnt++ is performed, the stencil won't reference cellset[128][j]
                cnt=0;
            }
        }
    }
    else if(i-3 <0){ //Stencil near the I=0 border
        int cnt=i+3;
        int n=6;
        while(n<7 && n>-1){ // Added -1 term to break the loop
            stencil[n--]=cellset[cnt--][j]; //I'm not sure n-- accounts for stencil[0] or at i=0 since stencil[2] = cellset[128][j] => nan Since cell 128 is a ghost cell
            if(cnt<0){ 
                cnt=cellset.size()-2; // ************** Changed this from cnt>=cellset.size()-1 so the stencil won't reference cellset[128][j], it will reference cellset[127]
            }
        }
    }
    return stencil;
}

vector<cellState> stencilNS(grid & grd, vector< vector<cellState> > & cellset, int i, int j){
    vector<cellState> stencil(7);
    
	if(j+3<cellset[0].size() && j-3 >=0){//Interior stencil 
        for(int n=0; n<7; n++){
            stencil[n]=cellset[i][j+n-3];
        }
    }
    else if(j+3>=cellset[0].size()){ // ******** CHANGED to j+3 >= size() ------------- j+3<cellset[0].size() was invalid since this would skip over 61. 61+3 = 64>64
        // INLET / OUTLET Stencil
        int cnt=j-3;
        int n=0;
        while(n<7){
            if(cnt>=grd.M-1){
                stencil[n++]=cellState(); //invalid cells across the boundary
            }
            else{
                stencil[n++]=cellset[i][cnt++];
            }
        }
    }
    else if(j-3 < 0){ //Stencil near the airfoil
        // AIRFOIL stencil
        int cnt=j+3; // Changed this from i+3 to j+3 from pattern matching. Is this right?
        bool turnaround=false;
        int n=6;
        while(n<7 && n>-1){ // Added -1 term to break the loop
            if(!turnaround){ // if turnaround is true
                stencil[n--]=cellset[i][cnt--]; // stencil[6] = cellset[i][3] -> stencil[5] = cellset[i][2] -> ... -> stencil[3] = cellset[i][0]
                if(cnt<0){
                    turnaround=true;
                    cnt=0;
                }
			}
			else if (n < 0) {

				std::cout << "n < 0 in StencilNS";
				system("pause");

			} else {
                stencil[n--]=cellset[i][cnt++]; //stencil[2] = cellset[i][0] -> stencil[1] = cellset[i][1] -> stencil[0] = cellset[i][2] --------------- CLEVER!!!!
            }
        }
    }
    return stencil;
}
