#include<iostream>
#include<fstream>
#include<sstream>
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
string solutionString(int i, int j);

int main(){
    cout << "Loading Initial Conditions\n";
    loadICs();
    writeGrid();

	//cout << "Please make the console full screen"<<endl;
	//system("pause");
    int writemodulo=1;//write solution every # of time steps.
    int tmax=5;
    for (int t=0 ; t<=tmax; t++){
		cout << "Time Step:\t"<<t<<"\n";
		//significantly speed up code by not writing every timestep (writing is slow)
		if(t==0 || t==tmax || t%writemodulo==0 || writemodulo <0){ //write solution set for t=0 , on write modulo steps, and in case of bad write modulo
            writeSolutionStep(t);
		}
        //RK4(grd, cellset, CFL);
        singleStepIntegration(grd, cellset, CFL);
    }
}

void loadICs(){
    cellset.resize(grd.N-1, std::vector<cellState>(grd.M-1));
    for(int i=0; i<grd.N-1;i++){
        for(int j=0; j<grd.M-1;j++){
            cellset[i][j].redefine(rho_ref,rho_ref*M_ref*c_ref,0,rhoE_ref,gamma,cv,i,j);
            if(i==80 && j>10 && j< 30){ // ************************* Fix me later **************************
                cellset[i][j].rhois(.5*rho_ref);
            }
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
			fout << grd.xCorner[i][j] << ' ' << grd.yCorner[i][j] << endl;
		}
	}

	fout.close();

	fout.open("GridCenters.dat");
	fout << "VARIABLES = \"Xctr\", \"Yctr\", \"xSds\", \"ySds\", \"xSnorm\",\"ySnorm\",";
	fout <<" \"xWds\", \"yWds\",\"xWnorm\", \"yWnorm\", \"xInorm\", \"yInorm\", \"xJnorm\", \"yJnorm\" \n";
	fout << "ZONE T = \"Grid center vals" << " \", I = " << grd.N << " , J = " << grd.M - 1 << ", F=POINT \n\n";
	for (int j = 0; j < grd.M - 1; j++) {
		for (int i = 0; i < grd.N - 1; i++) {
			fout << grd.xCenter[i][j] << ' ' << grd.yCenter[i][j] << ' ' << grd.xSdeltas[i][j] << ' ' << grd.ySdeltas[i][j] << ' ' << grd.xSnorm[i][j] << ' ';
			fout << grd.ySnorm[i][j] << ' ' << grd.xWdeltas[i][j] << ' ' << grd.yWdeltas[i][j] << ' ' << grd.xWnorm[i][j] << ' ';
			fout << grd.yWnorm[i][j] << ' ' << grd.xInorm[i][j] << ' ' << grd.yInorm[i][j] << ' ' << grd.xJnorm[i][j] << ' ' << grd.yJnorm[i][j] << endl;
		}
		fout << grd.xCenter[0][j] << ' ' << grd.yCenter[0][j] << ' ' << grd.xSdeltas[0][j] << ' ' << grd.ySdeltas[0][j] << ' ' << grd.xSnorm[0][j] << ' ';
		fout << grd.ySnorm[0][j] << ' ' << grd.xWdeltas[0][j] << ' ' << grd.yWdeltas[0][j] << ' ' << grd.xWnorm[0][j] << ' ';
		fout << grd.yWnorm[0][j] << ' ' << grd.xInorm[0][j] << ' ' << grd.yInorm[0][j] << ' ' << grd.xJnorm[0][j] << ' ' << grd.yJnorm[0][j] << endl;
	}

	fout.close();

	fout.open("EdgesAtCorners.dat");
	fout << "VARIABLES = \"xCorner\", \"yCorner\", \"xSds\", \"ySds\", \"xSnorm\", \"ySnorm\" \n";
	fout << "ZONE T = \"Outer edge" << " \", I = " << grd.N << " , J = " << grd.M << ", F=POINT \n\n";
	for (int j = 0; j < grd.M; j++) {
		for (int i = 0; i < grd.N - 1; i++) {
			//fout << "# i = " << i+1 << ", j = " << j+1 << endl;
			fout << grd.xCorner[i][j] << ' ' << grd.yCorner[i][j] << ' ' << grd.xSdeltas[i][j] << ' ' << grd.ySdeltas[i][j] << ' ' << grd.xSnorm[i][j] << ' ';
			fout << grd.ySnorm[i][j] << endl;
		}

		fout << grd.xCorner[0][j] << ' ' << grd.yCorner[0][j] << ' ' << grd.xSdeltas[0][j] << ' ' << grd.ySdeltas[0][j] << ' ' << grd.xSnorm[0][j] << ' ';
		fout << grd.ySnorm[0][j] << endl;

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
	fout.open("SolutionFile.dat", mode);
	if(!fout.good()){
        cout << "ERROR OPENING SOLUTION FILE\n";
	}
	fout << "VARIABLES = \"X\", \"Y\", \"Speed\", \"P\", \"M\",\"H\", \"S\", \"Xvel\", \"Yvel\",";
	fout << "\"F0\", \"F1\", \"F2\", \"F3\", \"G0\", \"G1\", \"G2\", \"G3\", \"rho\", ";
	fout << "\"RES0\", \"RES1\", \"RES2\", \"RES3\" \n";
	fout << "ZONE T = \"CELL CENTERS AT TIMESTEP " << t << " \", I = " << grd.N << " , J = " << grd.M - 1 << ", F=POINT \n\n";

	for (int j = 0; j < grd.M - 1; j++) {
		for (int i = 0; i < grd.N - 1; i++) {
            fout << solutionString(i,j);
		}
        fout << solutionString(0,j);
	}
//		for (int j = 0; j < grd.M - 1; j++) {
//            fout << grd.xCenter[0][j] << ' ' << grd.yCenter[0][j] << ' ' << cellset[0][j].speed() << ' ' << cellset[0][j].P() << ' ' << cellset[0][j].M() << ' ';
//            fout << cellset[0][j].H() << ' ' << cellset[0][j].S() << ' ' << cellset[0][j].U() << std::endl;
//		}
    fout << endl << "TEXT X = " << grd.N - 1 << ", Y = " << grd.M - 1 << ", T = \"Timestep = " << t << " \", F = COURIER, CS = FRAME, H = 2, ZN = " << t << endl << endl;
    //system("pause");
    fout.close();
}

string solutionString(int i, int j){
    stringstream fout;
	vector<double> residual(4,0.0);

    fout << endl << "# i = " << i << ", j = " << j << endl << endl;

    fout << grd.xCenter[i][j] << ' ' << grd.yCenter[i][j] << ' ' << cellset[i][j].speed() << ' ' << cellset[i][j].P() << ' ' << cellset[i][j].M() << ' ';
    fout << cellset[i][j].H() << ' ' << cellset[i][j].S() << ' ' << cellset[i][j].U() << ' ' << cellset[i][j].V() << ' ';
    fout << cellset[i][j].F1() << ' ' << cellset[i][j].F2() << ' ' << cellset[i][j].F3() << ' ' << cellset[i][j].F4() << ' ';
    fout << cellset[i][j].G1() << ' ' << cellset[i][j].G2() << ' ' << cellset[i][j].G3() << ' ' << cellset[i][j].G4() << ' ';
    fout << cellset[i][j].rho() << ' ';
    residual = Residuals(grd, cellset, i, j);
    fout << residual[0] << ' ' << residual[1] << ' ' << residual[2] << ' ' << residual[3] << endl;
    return fout.str();
}
