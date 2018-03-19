
//	Solver.h
//
//	This class can be used for solving systems of linear equation.
//  ElMtxL and BElMtxL are the lists of internal (triangular)
//	and boundary element matrices that define the system to be
//	solved.  These lists must be generated outside this class
//	and then passed to the solver of choice.
//	There are currently three linear solvers in this class:
//	1) a conjugate gradients solver
//	2) a preconditioned conjugate gradients solver
//	3) a preconditioned GMRES (Generalized Minimal RESisdual) solver
//
//	nVars is the number of variables that are being solved at each node
//
//	For the GRMES solver...
//	m = No. of Krylov Vectors or No. of steps before restart 
//	k = No. of gmres iterations
//	rec = convergence tolerance
////////////////////////////////////////////////////////

#ifndef SOLVER_H
#define SOLVER_H

#include "HabitatTIN.h"

class Solver
{
	protected:
		HabitatTIN *theMeshP;
		double Pen;

		void matrixVectorMult(double *InVec, double *OutVec, ItemList *ElMtxL, ItemList *BElMtxL, int nVars);
		void assembleRHS(double * RHS, ItemList *ElMtxL, ItemList *BElMtxL, int nVars);
		void MInverse(double *M, ItemList *ElMtxL, ItemList *BElMtxL, int nVars);
		double calcPenalty();	
	
	public:
		Solver(HabitatTIN* meshP);
		void preconditionedCGs(ItemList *ElMtxL, ItemList *BElMtxL, int nVars);
		void conjugateGradients(ItemList *ElMtxL, ItemList *BElMtxL, int nVars);
		void preconditionedGMRES(ItemList *ElMtxL, ItemList *BElMtxL, int nVars, int m, int k, double rec);

};

#endif
