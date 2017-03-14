#ifndef EULERRELATIONS_H_INCLUDED
#define EULERRELATIONS_H_INCLUDED
#include "cellState.h"

//X density flux contributions from a boundary value. Takes two cells and should make an averaged cell at the boundary by calling the cellState(cellState, cellState) constructor
double F1(cellState cellA, cellState cellB){
    cellState temp=cellState(cellA, cellB);
    double dL = cellA.eastXnorm()*cellA.eastYlength();
    return temp.rho()*temp.U()*dL;
}
double F2(cellState cell, cellState cellB);
double F3(cellState cell, cellState cellB);
double F4(cellState cell, cellState cellB);

//Y density flux contributions from a boundary value. Analogous to above
double G1(cellState cell, cellState cellB);
double G2(cellState cell, cellState cellB);
double G3(cellState cell, cellState cellB);
double G4(cellState cell, cellState cellB);

#endif // EULERRELATIONS_H_INCLUDED
