#ifndef CELLSTATE_H
#define CELLSTATE_H
#include <math.h>
#include <grid.h>

class cellState{
    //Cell State hold conserved properties of a cell as well as derived and thermo properties.
    //Variables xi, and eta link the cell to an independant grid object
    //Maybe this should include a pointer to the grid this cell is associated with??
    //That way functions the cell can fetch grid parameters, without needing to hold redundant data.
    private:
        double cell_rho, cell_rhoU, cell_rhoV, cell_rhoE, cell_gamma, cell_cv; //cell state properties
        int xi, eta; //(xi, eta) location of cell center;

    public:
        //constructors in cellState.cpp
        cellState(){
            cell_rho=0;
            cell_rhoU=0;
            cell_rhoV=0;
            cell_rhoE=0;
            cell_gamma=0;
            cell_cv=0;
        }
        cellState(double a, double b, double c, double d, double e, double f, int i, int j){ //constructs cell with a grid context
            cell_rho=a;
            cell_rhoU=b;
            cell_rhoV=c;
            cell_rhoE=d;
            cell_gamma=e;
            cell_cv=f;
            xi=i;
            eta=j;
        }
        cellState(double a, double b, double c, double d, double e, double f){ //makes an anonymous cell state manually
            cell_rho=a;
            cell_rhoU=b;
            cell_rhoV=c;
            cell_rhoE=d;
            cell_gamma=e;
            cell_cv=f;
        }
        cellState(cellState a, cellState b){//makes an anonymous cell state from the average of two cells
            cell_rho=(a.rho()+b.rho())/2;
            cell_rhoU=(a.rhoU()+b.rhoU())/2;
            cell_rhoV=(a.rhoV()+b.rhoV())/2;
            cell_rhoE=(a.rhoE()+b.rhoE())/2;
            cell_gamma=(a.gamma()+b.gamma())/2;
            cell_cv=(a.cv()+b.cv())/2;
        }
        //conserved quantities and thermo properties
        double rho() {return cell_rho;}
        double rhoU() {return cell_rhoU;}
        double rhoV() {return cell_rhoV;}
        double rhoE() {return cell_rhoE;}
        double gamma() {return cell_gamma;}
        double cv() {return cell_cv;}
        //derived properties
        double P(){
            return (cell_rhoE-.5*(pow(cell_rhoU,2)+pow(cell_rhoV,2))/cell_rho)*(cell_gamma-1);
        }
        double C() {return sqrt(cell_gamma*P()/cell_rho);}
        double H() {return (cell_rhoE+P())/cell_rho;}
        double T() {return P()/cell_rho/cell_cv;}
        double S() {return cell_cv*(log(T())-(cell_gamma-1)*log(cell_rho));}
        double U() {return cell_rhoU/cell_rho;}
        double V() {return cell_rhoV/cell_rho;}
        double speed() {return sqrt(pow(U(),2)+pow(V(),2));}
        double M() {return speed()/C();}
        double F1() {return cell_rhoU;}
        double F2() {return cell_rhoU * U() + P();}
        double F3() {return cell_rhoU * V();}
        double F4() {return cell_rhoU * H();}
        double G1() {return cell_rhoV;}
        double G2() {return cell_rhoV * V() + P();}
        double G3() {return cell_rhoV * U();}
        double G4() {return cell_rhoV * H();}
        //variable manipulation functions returns the value which has been set (Maybe they should return the THIS or void)
        double rhois(double a){
            cell_rho=a;
            return cell_rho;
        }
        double rhoUis(double a){
            cell_rhoU=a;
            return cell_rhoU;
        }
        double rhoVis(double a){
            cell_rhoV=a;
            return cell_rhoV;
        }
        double rhoEis(double a){
            cell_rhoE=a;
            return cell_rhoE;
        }
        void redefine(double a, double b, double c, double d, double e, double f, int i, int j){ //redefines cell with a grid context
            cell_rho=a;
            cell_rhoU=b;
            cell_rhoV=c;
            cell_rhoE=d;
            cell_gamma=e;
            cell_cv=f;
            xi=i;
            eta=j;
        }
        void redefine(double a, double b, double c, double d, double e, double f){
            cell_rho=a;
            cell_rhoU=b;
            cell_rhoV=c;
            cell_rhoE=d;
            cell_gamma=e;
            cell_cv=f;
        }
        void redefine(double a, double b, double c, double d){
            cell_rho=a;
            cell_rhoU=b;
            cell_rhoV=c;
            cell_rhoE=d;
        }
        //Following this guide on assignment operators overloading: https://www.tutorialspoint.com/cplusplus/assignment_operators_overloading.htm
        void operator = (const cellState &C ) { //Verified to work - Mike 3_10
            cell_rho=C.cell_rho;
            cell_rhoU=C.cell_rhoU;
            cell_rhoV=C.cell_rhoV;
            cell_rhoE=C.cell_rhoE;
            cell_gamma=C.cell_gamma;
            cell_cv=C.cell_cv;
            xi=C.xi;
            eta=C.eta;
        }
        int i(){return xi;}
        int j(){return eta;}
        //grid based functions
};
#endif // CELLSTATE_H
