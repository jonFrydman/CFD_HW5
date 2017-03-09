//#include "cellState.h"
//
//cellState::cellState(double a, double b, double c, double d, double e, double f, int i, int j){//this makes a specific cell in the mesh
//    cell_rho=a;
//    cell_rhoU=b;
//    cell_rhoV=c;
//    cell_rhoE=d;
//    cell_gamma=e;
//    cell_cv=f;
//    xi=i;
//    eta=j;
//}
//cellState::cellState(double a, double b, double c, double d, double e, double f){ //makes an anonymous cell state manually
//    cell_rho=a;
//    cell_rhoU=b;
//    cell_rhoV=c;
//    cell_rhoE=d;
//    cell_gamma=e;
//    cell_cv=f;
//}
//cellState::cellState(cellState a, cellState b){//makes an anonymous cell state from the average of two cells
//    cell_rho=(a.rho()+b.rho())/2;
//    cell_rhoU=(a.rhoU()+b.rhoU())/2;
//    cell_rhoV=(a.rhoV()+b.rhoV())/2;
//    cell_rhoE=(a.rhoE()+b.rhoE())/2;
//    cell_gamma=(a.gamma()+b.gamma())/2;
//    cell_cv=(a.cv()+b.cv())/2;
//}
//
////cellState::~cellState(){
////    //DESTRUCTOR, Use for deleting instances of 'new'
////}
