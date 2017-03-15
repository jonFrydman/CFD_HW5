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

int main(){
//    cout <<rho_ref <<"\n";
//    cellState cell(rho_ref, 10, 0, rhoE_ref, gamma, cv);
//    cout<<cell.rhoU()<<"\t"<<P_ref<<"\n";
////    cell.rhou()=5;
//    cout <<cell.rhoUis(5);
//    cout << "Grid Size is: " << grd.N <<" x "<<grd.M<<"\t (N x M)\n";
//    for(int i=0;i<grd.N-1; i++){
//        cout << "West Norm:\t<" << grd.xWnorm[i][0]<<", "<<grd.yWnorm[i][0]<<">\t\tSouth Norm:\t<"<< grd.xSnorm[i][0]<<", "<<grd.ySnorm[i][0]<<">";
//        cout<<"\t\tW-magnitude:\t"<<sqrt(pow(grd.xWnorm[i][0],2)+pow(grd.yWnorm[i][0],2));
//        cout<<"\t\tS-magnitude:\t"<<sqrt(pow(grd.xSnorm[i][0],2)+pow(grd.ySnorm[i][0],2));
//        cout<<"\t\tCell Area:\t"<<grd.area[i][0]<<"\t\tAspect Ratio:\t"<<grd.xSside[i][0]/grd.xWside[i][0]<<"\t\tOrthogonality:\t"<<cos(grd.xWside[i][0]/grd.xSside[i][0])*180/pi<<" deg\n\n";
//    }
    loadICs();

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

//    cout << "Cell Velocity:\t<"<<cellset[10][10].U()<<", "<<cellset[10][10].V()<<">\t\tPressure:\t"<<cellset[10][10].P()<<"\n";
//    cout << "Engergy:\t"<<cellset[10][10].rhoE()<<"\n";
//    cout << "P_ref:\t"<<P_ref<<"\t\trhoE_ref:\t"<<rhoE_ref<<"\t\trhoU_ref:\t"<<rho_ref*speed_ref<<"\n";
//
//    cout << "Expected %P rise at front of airfoil:\t"<<100*abs(P_ref-cellState(rho_ref,0,0,rhoE_ref,gamma,cv).P())/P_ref<<"%\n";
//
//    vector<double> test = EastFlux(grd, cellset, 10,10);
//    cout<<"F*(AV)1:\t" << test[0] << "\t\tF*(AV)2:\t" << test[1] << "\t\tF*(AV)3:\t" << test[2] << "\t\tF*(AV)4:\t" << test[3];

	cout << "Please make the console full screen"<<endl;
	system("pause");
	
	fout.open("SolutionFile.dat");
	//fout << "FILETYPE = SOLUTION \n";


    for (int t=0 ; t<3; t++){

        cout << "Time Step:\t"<< t << "\t";
        cout << "Cell Velocity: <"<<cellset[10][10].U() <<", "<<cellset[10][10].V()<<">\tPressure: "<<cellset[10][10].P();
        cout<<"\tCell Tau: "<< Tau(grd,cellset[10][10],CFL)<<"\t";
        cout<<"\t" << "Energy: "<<cellset[10][10].rhoE()<<"\t" << "Entropy: "<< cellset[10][10].S() << "\t" << "Temp: " << cellset[10][10].T() << "\t" << "Enthalpy: " << cellset[10][10].H()<<"\t"
        <<"Expected %P rise at front of airfoil: "<<100*abs(P_ref-cellState(rho_ref,0,0,rhoE_ref,gamma,cv).P())/P_ref<<"%\n";
        cout << endl;



		fout << "VARIABLES = \"X\", \"Y\", \"Speed\", \"P\", \"M\",\"H\", \"S\", \"Xvel\" \n";
		fout << "ZONE T = \"CELL CENTERS AT TIMESTEP " << t << " \", I = " << grd.N - 1 << " , J = " << grd.M - 1 << ", F=POINT \n\n";
		for (int j = 0; j < grd.M - 1; j++) {
			for (int i = 0; i < grd.N - 1; i++) {

				fout << grd.xCenter[i][j] << ' ' << grd.yCenter[i][j] << ' ' << cellset[i][j].speed() << ' ' << cellset[i][j].P() << ' ' << cellset[i][j].M() << ' ' << cellset[i][j].H()
					<< ' ' << cellset[i][j].S() << ' ' << cellset[i][j].U() << std::endl;

			}
		}

		fout << endl << "TEXT X = " << grd.N - 1 << ", Y = " << grd.M - 1 << ", T = \"Timestep = " << t << " \", F = COURIER, CS = FRAME, H = 2, ZN = " << t << endl << endl;


		system("pause");
        cellset = RK4(grd, cellset, CFL);

    }

	fout.close();

}

void loadICs(){
    cellset.resize(grd.N
		, std::vector<cellState>(grd.M-1));
    for(int i=0; i<grd.N-1;i++){
        for(int j=0; j<grd.M-1;j++){
            cellset[i][j].redefine(rho_ref,rho_ref*M_ref*c_ref,0,rhoE_ref,gamma,cv,i,j);
        }
    }
}
