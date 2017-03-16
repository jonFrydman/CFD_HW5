#include<iostream>
#include<fstream>
#include "referenceParameters.h"
#include "integration.h"
#include "eulerFlux.h"
#include "cellState.h"
#include "grid.h"
#include "numericalFlux.h"

using namespace std;

const char XFILENAME[]="Grid_X_Points.csv";
const char YFILENAME[]="Grid_Y_Points.csv";
const double pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062;

grid grd=grid(XFILENAME, YFILENAME);
vector< vector<cellState> > cellset;
void loadICs();
void writeGrid();
void writeSolutionStep(int t);

int main(){
    cout << "Loading Initial Conditions\n";
    loadICs();
    writeGrid();

	//cout << "Please make the console full screen"<<endl;
	//system("pause");

    for (int t=0 ; t<2; t++){
        writeSolutionStep(t);
        cout << "Time Step:\t"<<t<<"\n";
        cellset = RK4(grd, cellset, CFL);
    }
}

void loadICs(){
    cellset.resize(grd.N-1, std::vector<cellState>(grd.M-1));
    for(int i=0; i<grd.N-1;i++){
        for(int j=0; j<grd.M-1;j++){
            cellset[i][j].redefine(rho_ref,rho_ref*M_ref*c_ref,0,rhoE_ref,gamma,cv,i,j);
        }
    }
}
void writeGrid(){
	ofstream fout;
	fout.open("GridFile.dat");
	fout << "TITLE = NACA GRID \n";
	fout << "FILETYPE = GRID \n";
	fout << "VARIABLES = \"X\", \"Y\" \n";
	fout << "ZONE I=" << grd.N << " , J=" << grd.M << ", F=POINT \n";
	for (int j = 0; j < grd.M; j++) {
		for (int i = 0; i < grd.N; i++) {
			fout << grd.xCorner[i][j] << ' ' << grd.yCorner[i][j] << std::endl;
		}
	}
	fout.close();
}

void writeSolutionStep(int t){
    ofstream fout;
    ios_base::openmode mode;
    if(t==0){//if this is the first time step create a new file
        mode =ofstream::trunc;
    }
    else{//if this is not the first time step just append to the existing file
        mode=ofstream::app;
    }
	fout.open("SolutionFile.d   at", mode);
	if(!fout.good()){
        cout << "ERROR OPENING SOLUTION FILE\n";
	}

    fout << "VARIABLES = \"X\", \"Y\", \"Speed\", \"P\", \"M\",\"H\", \"S\", \"Xvel\", \"Yvel\", \"F0\", \"G0\", \"rho\", \"nSx\", \"nSy\", \"nWx\", \"nWy\", \n";
    fout << "ZONE T = \"CELL CENTERS AT TIMESTEP " << t << " \", I = " << grd.N - 1 << " , J = " << grd.M - 1 << ", F=POINT \n\n";
    for (int j = 0; j < grd.M - 1; j++) {
        for (int i = 0; i < grd.N-1; i++) {
            fout << grd.xCenter[i][j] << ' ' << grd.yCenter[i][j] << ' ' << cellset[i][j].speed() << ' ' << cellset[i][j].P() << ' ' << cellset[i][j].M() << ' ';
            fout << cellset[i][j].H() << ' ' << cellset[i][j].S() << ' ' << cellset[i][j].U() << ' ' << cellset[i][j].V() << ' ';
            fout << cellset[i][j].F1() << ' ' << cellset[i][j].G1() << ' ' << cellset[i][j].rho() << ' ';
            fout << grd.xSnorm[i][j] << ' ' << grd.ySnorm[i][j] << ' ' << grd.xWnorm[i][j] << ' ' << grd.yWnorm[i][j] << endl;
        }
    }
//		for (int j = 0; j < grd.M - 1; j++) {
//            fout << grd.xCenter[0][j] << ' ' << grd.yCenter[0][j] << ' ' << cellset[0][j].speed() << ' ' << cellset[0][j].P() << ' ' << cellset[0][j].M() << ' ';
//            fout << cellset[0][j].H() << ' ' << cellset[0][j].S() << ' ' << cellset[0][j].U() << std::endl;
//		}
    fout << endl << "TEXT X = " << grd.N - 1 << ", Y = " << grd.M - 1 << ", T = \"Timestep = " << t << " \", F = COURIER, CS = FRAME, H = 2, ZN = " << t << endl << endl;
    //system("pause");
    fout.close();
}
