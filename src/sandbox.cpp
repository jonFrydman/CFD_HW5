#include<iostream>
#include<fstream>
#include "referenceParameters.h"
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

 for (int t;t<10;t++){

    cout << "Time Step:\t"<< t << endl;
    cout << "Cell Velocity:\t<"<<cellset[10][10].U()<<", "<<cellset[10][10].V()<<">\t\tPressure:\t"<<cellset[10][10].P()<<"\n";
    cout << "Energy:\t"<<cellset[10][10].rhoE()<<"\n";
    cout << "P_ref:\t"<<P_ref<<"\t\trhoE_ref:\t"<<rhoE_ref<<"\t\trhoU_ref:\t"<<rho_ref*speed_ref<<"\n";

    cout << "Expected %P rise at front of airfoil:\t"<<100*abs(P_ref-cellState(rho_ref,0,0,rhoE_ref,gamma,cv).P())/P_ref<<"%\n";

    vector<double> test = EastFlux(grd, cellset, 10,10);
    cout<<"\t\tF*(AV)1:\t" << test[0] << endl << "\t\tF*(AV)2:\t" << test[1]<< endl  << "\t\tF*(AV)3:\t" << test[2]<< endl << "\t\tF*(AV)4:\t" << test[3]<< endl;

    system("pause");


}

return 0;
}

void loadICs(){
    cellset.resize(grd.N-1, std::vector<cellState>(grd.M-1));
    for(int i=0; i<grd.N-1;i++){
        for(int j=0; j<grd.M-1;j++){
            cellset[i][j].redefine(rho_ref,rho_ref*M_ref*c_ref,0,rhoE_ref,gamma,cv,i,j);
        }
    }
}



