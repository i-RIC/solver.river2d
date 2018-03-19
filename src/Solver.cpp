
//		Solver.cpp

//		Solver class function definitions

//#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <ctype.h>
#include "Solver.h"
#include "ElementMtx.h"
#include "Element.h"
#include "Shallow.h"
using namespace std;

Solver::Solver(HabitatTIN *meshP)
{
	theMeshP = meshP;
	Pen = 0;
}

void Solver::matrixVectorMult(double *InVec, double *OutVec, ItemList *ElMtxL, ItemList *BElMtxL, int nVars)
{

	int n, nn, inVectrow, outVectrow;
	ElementMtx *elMP, *BelMP;
	Element *elP, *belP;
	Node *inVectnP, *outVectnP;
	int i, j , ii, jj;

	n = theMeshP->Nnodes() * nVars;

	//initialize the out vector

	for (i = 0; i < n; i ++)
	{
		OutVec[i] = 0;
	}

	elMP = (ElementMtx*)ElMtxL->firstItem();
	BelMP = (ElementMtx*)BElMtxL->firstItem();

	//interior element contributions

	while (elMP != NULL)
	{
		elP = elMP->getElement();
		nn = elP->numNodes();
		for (i = 0; i < nn; i++)
		{
			outVectnP = elP->getNode(i);
			outVectrow = outVectnP->getIndex()*nVars;
			for (j = 0; j < nn;  j++)
			{
				inVectnP = elP->getNode(j);
				inVectrow = inVectnP->getIndex()*nVars;
				for( ii = 0; ii < nVars; ii++)
				{
					for(jj = 0; jj < nVars; jj++)
					{
						OutVec[outVectrow + ii] += elMP->KE(i*nVars+ii, j*nVars+jj) * InVec[inVectrow + jj];
					}
				}
			}
		}
	elMP = (ElementMtx*)ElMtxL->nextItem();
	}

	//boundary element contributions

	while (BelMP != NULL)
	{
		belP = BelMP->getElement();
		nn = belP->numNodes();
		for (i = 0; i < nn; i++)
		{
			outVectnP = belP->getNode(i);
			outVectrow = outVectnP->getIndex()*nVars;
			for (j = 0; j < nn;  j++)
			{
				inVectnP = belP->getNode(j);
				inVectrow = inVectnP->getIndex()*nVars;
				for( ii = 0; ii < nVars; ii++)
				{
					for(jj = 0; jj < nVars; jj++)
					{
						OutVec[outVectrow + ii] += BelMP->KE(i*nVars+ii, j*nVars+jj) * InVec[inVectrow + jj] * Pen;
					}
				}
			}
		}
	BelMP = (ElementMtx*)BElMtxL->nextItem();
	}
}

void Solver::assembleRHS(double *RHS, ItemList *ElMtxL, ItemList *BElMtxL, int nVars)
{
	int n, nn, outVectrow;
	ElementMtx *elMP, *BelMP;
	Element *elP, *belP;
	Node * outVectnP;
	int i, ii;

	n = theMeshP->Nnodes() * nVars;

	//initialize the out vector

	elMP = (ElementMtx*)ElMtxL->firstItem();
	BelMP = (ElementMtx*)BElMtxL->firstItem();
	
	//initialize the RHS vector

	for (i = 0; i < n; i ++) RHS[i] = 0;

	//interior element contributions

	while (elMP != NULL)
	{
		elP = elMP->getElement();
		nn = elP->numNodes();
		for (i = 0; i < nn; i++)
		{
			outVectnP = elP->getNode(i);
			outVectrow = outVectnP->getIndex()*nVars;
			for(ii = 0; ii < nVars; ii++)
			{
				RHS[outVectrow + ii] += elMP->FE(i*nVars+ii);
			}
		}
		elMP = (ElementMtx*)ElMtxL->nextItem();
	}

	//boundary element contributions

	while (BelMP != NULL)
	{		
		belP = BelMP->getElement();
		nn = belP->numNodes();
		for(i = 0; i < nn; i++)
		{
			outVectnP = belP->getNode(i);
			outVectrow = outVectnP->getIndex()*nVars;
			for(ii = 0; ii < nVars; ii++)
			{
				RHS[outVectrow + ii] += BelMP->FE(i*nVars+ii)*Pen;
			}
		}
		BelMP = (ElementMtx*)BElMtxL->nextItem();
	}
}

void Solver::MInverse(double *PreCon, ItemList *ElMtxL, ItemList *BElMtxL, int nVars)  
{
	int n, nn, outVectrow;
	ElementMtx *elMP, *BelMP;
	Element *elP, *belP;
	Node * outVectnP;
	int i, ii;

	n = theMeshP->Nnodes() * nVars;

	//initialize the out vector

	elMP = (ElementMtx*)ElMtxL->firstItem();
	BelMP = (ElementMtx*)BElMtxL->firstItem();
	
	//initialize the RHS vector

	for (i = 0; i < n; i ++) PreCon[i] = 0;

	//interior element contributions

	while (elMP != NULL)
	{
		elP = elMP->getElement();
		nn = elP->numNodes();
		for (i = 0; i < nn; i++)
		{
			outVectnP = elP->getNode(i);
			outVectrow = outVectnP->getIndex()*nVars;
			for(ii = 0; ii < nVars; ii++)
			{
				PreCon[outVectrow + ii] += elMP->KE(i*nVars+ii, i*nVars+ii);
			}
		}
		elMP = (ElementMtx*)ElMtxL->nextItem();
	}

	//boundary element contributions
	
	while (BelMP != NULL)
	{		
		belP = BelMP->getElement();
		nn = belP->numNodes();
		for(i = 0; i < nn; i++)
		{
			outVectnP = belP->getNode(i);
			outVectrow = outVectnP->getIndex()*nVars;
			for(ii = 0; ii < nVars; ii++)
			{
				PreCon[outVectrow + ii] += BelMP->KE(i*nVars+ii, i*nVars+ii)*Pen;
			}
		}
		BelMP = (ElementMtx*)BElMtxL->nextItem();
	}

	for (i = 0; i < n; i ++) PreCon[i] = 1.0/PreCon[i];
}

double Solver::calcPenalty()
{

	double penalty;
	double QIn = 0.0;
	RivBound *bP = (RivBound *)theMeshP->firstSeg();
	while(bP != NULL) 
	{
		if(bP->getBcCode() == 1)
		QIn += bP->getBcValue() * bP->length();
		bP = (RivBound *)theMeshP->nextSeg();
	}

	penalty = QIn * 100000;

	if (penalty < 1000000) penalty = 1000000;		// Added by CDA April 25, 2002 to ensure that the penalty has a large value even if the discarge is zero.

	return penalty;

}

void Solver::conjugateGradients(ItemList *ElMtxL, ItemList *BElMtxL, int nVars)
{
	// the conjugate gradients solution code

	double *X, *B, *Q, *D, *R;
	int i, imax, j, n;
	double deltaNew = 0, deltaOld = 0, delta0, err, err2delta0, alpha, beta, dTq;

	
	X = B = Q = D = R = NULL;

	n = theMeshP->Nnodes() * nVars;

	X = new double[n];  //had difficulty allocating these anywhere except here
	B = new double[n];
	Q = new double[n];
	D = new double[n];
	R = new double[n];

	//setting the initial conditions in the X vector

	theMeshP->initSolutionVector(X);

	//initializing variables for the iterative solution
	
	i = 0;
	imax = 2*n;		//theortically this method should converge in n iterations or less 
	err = 0.0000001;
	Pen = calcPenalty();
	matrixVectorMult(X,Q,ElMtxL, BElMtxL, nVars); // {Q} = [A]{X}
	assembleRHS(B,ElMtxL, BElMtxL, nVars);
	for (j = 0; j < n; j++) 
	{
		R[j] = B[j] - Q[j];
		D[j] = R[j];
		deltaNew += R[j]*R[j];
	}
	delta0 = deltaNew;
	err2delta0 = err*err*delta0;

	//the iterative solution

	while ( (i < imax) && (deltaNew > err2delta0) )
	{
		matrixVectorMult(D,Q,ElMtxL, BElMtxL, nVars);							// {Q} = [A]{D}
		dTq = 0;										// sum of {D[i]}*{Q[i]}		
		for(j = 0; j < n; j++) dTq = dTq + D[j]*Q[j];
		alpha = deltaNew / dTq;
		for(j = 0; j < n; j++) X[j] += alpha*D[j];		// {X} = {X} + alpha*{D}
		if ((i % 50) == 0)								// if divisable by 50
		{
			matrixVectorMult(X,Q,ElMtxL, BElMtxL, nVars);						// {Q} = [A]{X}
			for(j = 0; j < n; j++) R[j] = B[j] - Q[j];  // {R} = {B} - {Q}
		}
		else
			for(j = 0; j < n; j++) R[j] -= alpha*Q[j];  // {R} = {R} - alpha*{Q}
		deltaOld = deltaNew;
		deltaNew = 0;
		for(j = 0; j < n; j ++) deltaNew += R[j]*R[j];
		beta = deltaNew / deltaOld;
		for(j = 0; j < n; j ++) D[j] = R[j] + beta*D[j]; // {D} = {R} + beta*{D}
		i++;
	}

//	CString check;
//	check.Format("number of iterations = %d", i);
//	AfxMessageBox(check);

	//setting the nodal values of psi to the solution values

	theMeshP->updateSolToNodes(X);

	delete[] X;
	delete[] B;
	delete[] Q;
	delete[] D;
	delete[] R;

}

void Solver::preconditionedCGs(ItemList *ElMtxL, ItemList *BElMtxL, int nVars)
{

	double *X, *B, *Q, *D, *R, *Min, *S;
	int i, imax, j, n;
	double deltaNew = 0, deltaOld = 0, delta0, err, err2delta0, alpha, beta, dTq;

	X = B = Q = D = R = S = NULL;

	n = theMeshP->Nnodes() * nVars;

	X = new double[n];  //had difficulty allocating these anywhere except here
	B = new double[n];
	Q = new double[n];
	D = new double[n];
	R = new double[n];
	Min = new double[n];
	S = new double[n];

	//setting the initial conditions in the X vector

	theMeshP->initSolutionVector(X);

	//initializing variables for the iterative solution
	
	i = 0;
	imax = 2*n;	//theortically this method should converge in n iterations or less 
	err = 0.000000001;
	Pen = calcPenalty();
	matrixVectorMult(X, Q, ElMtxL, BElMtxL, nVars); // {Q} = [A]{X}
	assembleRHS(B, ElMtxL, BElMtxL, nVars);
	MInverse(Min, ElMtxL, BElMtxL, nVars);
	for (j = 0; j < n; j++) 
	{
		R[j] = B[j] - Q[j];
		D[j] = Min[j]*R[j];
		deltaNew += R[j]*D[j];
	}
	delta0 = deltaNew;
//	err2delta0 = err*err*delta0;
	err2delta0 = 0.0000001;

	//the iterative solution

	while ( (i < imax) && (deltaNew > err2delta0) )
	{
		matrixVectorMult(D, Q, ElMtxL, BElMtxL, nVars);		// {Q} = [A]{D}
		dTq = 0;											// sum of {D[i]}*{Q[i]}		
		for(j = 0; j < n; j++) dTq = dTq + D[j]*Q[j];
		alpha = deltaNew / dTq;
		for(j = 0; j < n; j++) X[j] += alpha*D[j];			// {X} = {X} + alpha*{D}
		if ((i % 50) == 0)									// if divisable by 50
		{
			matrixVectorMult(X, Q, ElMtxL, BElMtxL, nVars);	// {Q} = [A]{X}
			for(j = 0; j < n; j++) R[j] = B[j] - Q[j];		// {R} = {B} - {Q}
		}
		else
			for(j = 0; j < n; j++) R[j] -= alpha*Q[j];		// {R} = {R} - alpha*{Q}
		for(j = 0; j < n; j++) S[j] = Min[j]*R[j];			// {S} = [Min]{R}
		deltaOld = deltaNew;
		deltaNew = 0;
		for(j = 0; j < n; j ++) deltaNew += R[j]*S[j];
		beta = deltaNew / deltaOld;
		for(j = 0; j < n; j ++) D[j] = S[j] + beta*D[j];	 // {D} = {S} + beta*{D}
		i++;
	}

//	CString check;
//	check.Format("number of iterations = %d", i);
//	AfxMessageBox(check);

	//setting the nodal values to the solution values

	theMeshP->updateSolToNodes(X);

	delete[] X;
	delete[] B;
	delete[] Q;
	delete[] D;
	delete[] R;
	delete[] Min;
	delete[] S;

}

void Solver::preconditionedGMRES(ItemList *ElMtxL, ItemList *BElMtxL, int nVars, int m, int k, double rec)
{

	double *X, *B, *Q, *D, *R, *Min, *V, *H;
	double tem = 1, res = 0, res0, ccos, ssin, resCheck;
	int ii, it, i, j, n, nk, i0, im;
//	double err;
//	CString check;

	X = B = Q = D = R = H = NULL;

	n = theMeshP->Nnodes() * nVars; //number of unknowns

	X = new double[n];  //had difficulty allocating these anywhere except here
	B = new double[n];
	Q = new double[n];
	D = new double[n];
	R = new double[m+1];
	Min = new double[n];
	V = new double[n*m];
	H = new double[m*(m+1)];

	//setting the initial conditions in the X vector

	theMeshP->initSolutionVector(X);

	//initializing variables for the iterative solution

	Pen = calcPenalty();

	matrixVectorMult(X, Q, ElMtxL, BElMtxL, nVars); // {Q} = [A]{X}
	assembleRHS(B, ElMtxL, BElMtxL, nVars);
	MInverse(Min, ElMtxL, BElMtxL, nVars);

	for (ii = 0; ii < n; ii++) 
	{
		D[ii] = B[ii] - Q[ii];
		Q[ii] = Min[ii]*D[ii];
		res += Q[ii]*Q[ii];
	}

	res = sqrt(res);
	res0 = res;

	//the iterative solution

	for(it = 0; it < k; it++)
	{
		nk = m;
		
		//Arnoldi method (which uses the Gram-Schmidt orthogonalization method)
		for (ii = 0; ii < n; ii++) Q[ii] = Q[ii] / res;
		for (ii = 0; ii <= nk; ii++) R[ii] = 0;
		R[0] = res;

		for(j = 0; j < m; j ++)
		{
			for (ii = 0; ii < n; ii++) V[ii*m + j] = Q[ii];
			matrixVectorMult(Q, D, ElMtxL, BElMtxL, nVars);
			for (ii = 0; ii < n; ii++) Q[ii] = Min[ii]*D[ii];
			for(i = 0; i <= j; i++)
			{
				H[i*m + j] = 0;
				for(ii = 0; ii < n; ii++) H[i*m + j] += Q[ii]*V[ii*m + i];
				for(ii = 0; ii < n; ii++) Q[ii] -= H[i*m + j]*V[ii*m + i];
			}
			tem = 0;
			for (ii = 0; ii < n; ii++) tem += Q[ii]*Q[ii];
			tem = sqrt(tem);
			H[(j+1)*m + j] = tem;
			if (tem < rec)
			{
				nk = j + 1;
				goto l5;
			}
			for (ii = 0; ii < n; ii++) Q[ii] = Q[ii] / tem;
		}

		//triangularization
l5:		for (i = 0; i < nk; i++)
		{
			im = i + 1;
			tem = (1.0)/(sqrt(H[i*m + i]*H[i*m + i] + H[im*m + i]*H[im*m + i]));
			ccos = H[i*m + i]*tem;
			ssin = -H[im*m+i]*tem;
			for( j = i; j < nk; j++)
			{
				tem = H[i*m+j];
				H[i*m + j] = ccos * tem - ssin * H[im*m + j];
				H[im*m + j] = ssin * tem + ccos *H[im*m + j];
			}
			R[im] = ssin * R[i];
			R[i] = R[i]*ccos;
		}
		//solution of linear system
		for (i = (nk-1); i >= 0; i--)
		{
			R[i] = R[i]/H[i*m + i];
			for (i0 = (i - 1); i0 >= 0; i0--) R[i0] = R[i0] - H[i0*m + i]*R[i];
		}
		for (i = 0; i < nk; i++)
		{
			for(ii = 0; ii < n; ii++) 
				X[ii] += R[i]*V[ii*m + i];
		}
		//new residual and stopping tests
		matrixVectorMult(X,Q, ElMtxL, BElMtxL, nVars);
		res = 0;
		for(ii = 0; ii < n; ii++)
		{
			D[ii] = B[ii] - Q[ii];
			Q[ii] = Min[ii]*D[ii];
			res += Q[ii]*Q[ii];
		}
		res = sqrt(res);
		resCheck = res/res0;
//		check.Format("it = %d, res = %f", it, res);
//		AfxMessageBox(check);

		if(resCheck < rec) break;

	}


//	check.Format("number of iterations = %d, resCheck = %f res = %f, res0 = %f", it, resCheck, res, res0);
//	AfxMessageBox(check);

	//setting the nodal values to the solution values

	theMeshP->updateSolToNodes(X);


	delete[] X;
	delete[] B;
	delete[] Q;
	delete[] D;
	delete[] R;
	delete[] Min;
	delete[] V;
	delete[] H;

}
