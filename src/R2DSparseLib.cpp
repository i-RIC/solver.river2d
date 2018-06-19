//
//		R2DSparseLib.cpp
//
//		created:	Sept. 23, 2006
//		revised:	Oct. 1, 2007
//
//		Sparse Matrix methods for River2D programs
//

#include "R2DSparseLib.h"
#include "HabitatTIN.h"
#include "FE_PS.h"
//#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctype.h>
#include "Element.h"
#include "Shallow.h"

#define SE(A,B) ( *(ep.S + ep.n*(A) + (B)) )
#define FE(B)   ( *(ep.F + B) )
#define Vp(A,B) ( *(eqnsets[nset].Vp + eqnsets[nset].m*(A) + (B)) )
#define Hp(A,B) ( *(eqnsets[nset].Hp + eqnsets[nset].m*(A) + (B)) )

extern "C" struct	elmpointers		ep ;
extern "C" struct	pointers		gp;
extern "C" struct	control			N ;
extern "C" int		jacobianType;
extern "C" struct	eqnset			eqnsets[4] ;
extern "C" int		Nukns;
extern "C" double	uchange, maxChange;

extern "C" int		get_PGKeAnalJ(element *elmntp, element *theElp, int nvar, int code, int *ntf, int itnum);
extern "C" int		get_bKeAnalJ(belement *elmntp, belement *theElp, int nvar, int code, int *ntf, int itnum);
extern "C" int		get_PGKeNumJ(element *elmntp, element *theElp, int nvar, int code, int *ntf, int itnum);
extern "C" int		get_bKeNumJ(belement *elmntp, belement *theElp, int nvar, int code, int *ntf, int itnum);
extern "C" int		JacobiPC(double *Min, int itnum, int varnum, int elcode);

namespace R2DSparseLib {

	R2DSparseRow::R2DSparseRow(int row, int initialSize, int maxColNumber) {
		rowNumber = row;
		entries.reserve(initialSize);
		setZeroEntry(rowNumber);

/*		if (rowNumber > 1) 
			addValue(rowNumber-1, -1.0);
		if (rowNumber < maxColNumber)
			addValue(rowNumber+1, -1.0);
*/
	}

	R2DSparseRow::~R2DSparseRow() {}

	double R2DSparseRow::getValueByColumn(int column) {
		int foundFlag = 0;
		double foundValue;
		for(size_t i=0; i<entries.size(); i++){
			if(column == entries[i].columnIndex){
				foundValue = entries[i].value;
				foundFlag = 1;
			}
		}
		if(foundFlag == 0)
			foundValue = 0.0;
		return foundValue;
	}

	void R2DSparseRow::setEntry(int index, int ncolumn, double value){
		if(index < getRowSize()) {
			entries[index].columnIndex = ncolumn;
			entries[index].value = value;
		}
	}

	void R2DSparseRow::setValue(int ncolumn, double value){
		int newFlag = 1;
		for(size_t i=0;i<entries.size();i++){
			if(ncolumn == entries[i].columnIndex){
				entries[i].value = value;
				newFlag = 0;
				break;
			}
		}
		if(newFlag != 0){
			R2DSparseValue newEntry;
			newEntry.columnIndex = ncolumn;
			newEntry.value = value;
			entries.push_back(newEntry);
		}
	}

	double R2DSparseRow::getValue(int ncolumn){
		double value = 0.0;
		for(size_t i=0;i<entries.size();i++){
			if(ncolumn == entries[i].columnIndex){
				value = entries[i].value;
				break;
			}
		}
		return value;
	}

	void R2DSparseRow::setZeroEntry(int ncolumn){
		int newFlag = 1;
		for(size_t i=0;i<entries.size();i++){
			if(ncolumn == entries[i].columnIndex){
				entries[i].value = 0.0;
				newFlag = 0;
				break;
			}
		}
		if(newFlag != 0){
			R2DSparseValue newEntry;
			newEntry.columnIndex = ncolumn;
			newEntry.value = 0.0;
			entries.push_back(newEntry);
		}
	}

	void R2DSparseRow::addValue(int ncolumn, double value){
		int newFlag = 1;
		for(size_t i=0;i<entries.size();i++){
			if(ncolumn == entries[i].columnIndex){
				entries[i].value += value;
				newFlag = 0;
				break;
			}
		}
		if(newFlag != 0){
			R2DSparseValue newEntry;
			newEntry.columnIndex = ncolumn;
			newEntry.value = value;
			entries.push_back(newEntry);
		}
	}

	void R2DSparseRow::addRow(double factor, R2DSparseRow *pivotRow, int *indexArray) {

		R2DSparseValue fillEntry;

		for(int i=0;i<getRowSize();i++) {
			indexArray[entries[i].columnIndex] = i;
		}
		for(int j=0;j<pivotRow->getRowSize();j++) {
			int targetIndex = indexArray[pivotRow->entries[j].columnIndex];
			if( targetIndex > -1){
				entries[targetIndex].value += factor * pivotRow->entries[j].value;
			}
			else if(pivotRow->entries[j].value != 0){
				fillEntry.columnIndex = pivotRow->entries[j].columnIndex;
				fillEntry.value = pivotRow->entries[j].value;
				entries.push_back(fillEntry);
			}
		}
		for(int i=0;i<getRowSize();i++) {
			indexArray[entries[i].columnIndex] = -1;
		}
	}

	void R2DSparseRow::reduceRow(R2DSparseRow *pivotRow, int *indexArray) {
//		R2DSparseValue fillEntry;
		R2DSparseValue pivotEntry = pivotRow->entries[0];
		

		for(int i=0;i<getRowSize();i++) {
			indexArray[entries[i].columnIndex] = i;
		}
		if(indexArray[pivotEntry.columnIndex] > -1){
			int targetIndex = indexArray[pivotRow->entries[0].columnIndex];
			
			double factor = entries[targetIndex].value / pivotEntry.value;
			entries[targetIndex].value = factor;
			for(int j=1;j<pivotRow->getRowSize();j++) {
					if(pivotRow->entries[j].columnIndex > pivotRow->getRowNumber()) {
						
					int targetIndex = indexArray[pivotRow->entries[j].columnIndex];
					
					if( targetIndex > -1){
						entries[targetIndex].value -= factor * pivotRow->entries[j].value;	
					}	
				}
			}
		}
		for(int i=0;i<getRowSize();i++) {
			indexArray[entries[i].columnIndex] = -1;
		}
	}
	
	void R2DSparseRow::sweepRow(double *rhs){
		double sum = rhs[rowNumber-1];
		for(size_t i=1;i<entries.size();i++){
			if(entries[i].columnIndex < rowNumber){
				sum -= entries[i].value * rhs[entries[i].columnIndex-1];
			}
			
		}
		rhs[rowNumber-1] = sum; 
	}

	void R2DSparseRow::backSubRow(double *rhs){
		double sum = rhs[rowNumber-1];
		for(size_t i=1;i<entries.size();i++){
			if(entries[i].columnIndex > rowNumber){
				sum -= entries[i].value * rhs[entries[i].columnIndex-1];
			}
			
		}
		rhs[rowNumber-1] = sum / entries[0].value; 
	}

	double R2DSparseRow::vectorProduct(double *vector){
		double sum = 0.0;
		for(size_t i=0;i<entries.size();i++){
			sum += entries[i].value * vector[entries[i].columnIndex-1];
		}
		return sum; 
	}

	void R2DSparseRow::clearRow(){
		entries.clear();
	}


	R2DSparseMatrix::R2DSparseMatrix(int numRows, int rowSize) {
		nRows = numRows;
		rows = new R2DSparseRow*[nRows];
		indexArray = new int[nRows+1];
		for(int i=0; i<nRows; i++) {
			rows[i] = new R2DSparseRow(i+1,rowSize,nRows);
			indexArray[i+1] = -1;
		}
	}

	R2DSparseMatrix::~R2DSparseMatrix() {
		for(int i=0; i<nRows; i++) {
			delete rows[i];
		}
		delete [] rows;
		delete [] indexArray;
	}

	int R2DSparseMatrix::getNEntries() {
		int sum = 0;
		for(int i=0; i<nRows; i++) {
			sum += rows[i]->getRowSize();
		}
		return sum;
	}

	void R2DSparseMatrix::addValue(int row, int column, double value) {
		if(row<=nRows){
			rows[row-1]->addValue(column,value);
		}
	}

	void R2DSparseMatrix::setValue(int row, int column, double value) {
		if(row<=nRows){
			rows[row-1]->setValue(column,value);
		}
	}
	double R2DSparseMatrix::getValue(int row, int column){
		if(row<=nRows)
			return	rows[row-1]->getValue(column);
		else
			return 0.0;
	}

	void R2DSparseMatrix::addRow(int pivotRow, int targetRow, double factor){
		if (pivotRow <= nRows && targetRow <= nRows && pivotRow > 0 && targetRow >0) {
			R2DSparseRow *tRow = rows[targetRow-1];
			R2DSparseRow *pRow = rows[pivotRow-1];
			tRow->addRow(factor,pRow,indexArray);
		}
	}

	void R2DSparseMatrix::clearMatrix(){
		for(int i=0;i<nRows;i++){
			rows[i]->clearRow();
		}
	}

	void R2DSparseMatrix::addSubMatrix(int nr, int nc, const int *rownums, const int *colnums, double *sMatrix){
		R2DSparseRow newRow(rownums[0],nc,nRows);
		for(int i=0;i<nc;i++) {
			newRow.setZeroEntry(colnums[i]);
		}
		for(int i=0;i<nr;i++) {
			newRow.setRowNumber(rownums[i]);
			for(int j=0;j<nc;j++){
				newRow.setEntry(j,colnums[j],sMatrix [i*nc + j]);
			}
			rows[rownums[i]-1]->addRow(1.0,&newRow,indexArray);
		}
	}

	//void R2DSparseMatrix::addSubMatrix(ElementMatrix *sMatrix){
	//	int nr = sMatrix->getnRows();
	//	int nc = sMatrix->getnCols();
	//	int *rownums = sMatrix->getrowIndices();
	//	int *colnums = sMatrix->getcolIndices();
	//	R2DSparseRow newRow(rownums[0],nc,nRows);
	//	for(int i=0;i<nc;i++) {
	//		newRow.setZeroEntry(colnums[i]);
	//	}
	//	for(int i=0;i<nr;i++) {
	//		newRow.setRowNumber(rownums[i]);
	//		for(int j=0;j<nc;j++){
	//			newRow.setEntry(j,colnums[j],sMatrix->getMatrixValue(i,j));
	//		}
	//		rows[rownums[i]-1]->addRow(1.0,&newRow,indexArray);
	//	}
	//}

	void R2DSparseMatrix::reduceRow(int pivotRow, int targetRow){
		if (pivotRow <= nRows && targetRow <= nRows && pivotRow > 0 && targetRow >0) {
			R2DSparseRow *tRow = rows[targetRow-1];
			R2DSparseRow *pRow = rows[pivotRow-1];
			tRow->reduceRow(pRow,indexArray);
		}
	}

	void R2DSparseMatrix::LUDecompose(){
		R2DSparseRow *targetRow;
		R2DSparseRow *pivotRow;
		int iTarget;
		for(int iPivot=1; iPivot < nRows; iPivot++){
			pivotRow = rows[iPivot-1];
			for(int index=1; index<pivotRow->getRowSize(); index++) {
				iTarget = pivotRow->getColumnByIndex(index);
				if(iTarget > iPivot) {
					targetRow = rows[iTarget-1];
					targetRow->reduceRow(pivotRow,indexArray); 
				}
			}
		}
	}

	void R2DSparseMatrix::forwardSweep(double *rhs){
		for(int i=2; i <= nRows; i++){
			rows[i-1]->sweepRow(rhs);
		}
	}

	void R2DSparseMatrix::backSubstitute(double *rhs){
		for(int i=nRows; i >0; i--){
			rows[i-1]->backSubRow(rhs);
		}
	}

	void R2DSparseMatrix::vectorProduct(double *result, double *vector){
		for(int i=0; i <nRows; i++){
			result[i] = rows[i]->vectorProduct(vector);
		}
	}
	
	void R2DSparseMatrix::matrixVectorMultSparse(double *X, double *Q){
		int i,j;
		int nRows = getNRows();
		R2DSparseRow *currentRow;
		int currentRowSize;
		int currentColumnIndex, currentRowIndex;
		double currentValue;

		for(i=0;i<nRows;i++) Q[i]=0.0;

		for(i=0;i<nRows;i++)
		{
			currentRow = getRow(i+1);
			currentRowSize = currentRow->getRowSize();

			for (j=0;j<currentRowSize;j++)
			{
				currentColumnIndex = currentRow->getColumnByIndex(j);
				currentValue = currentRow->getValueByIndex(j);
				currentRowIndex = currentRow->getRowNumber();

				Q[currentRowIndex-1] += currentValue * X[currentColumnIndex-1];
			}
		}
	}
	
	R2DSparseMatrix::R2DSparseMatrix (const R2DSparseMatrix &K,int numRows) { // Copy constructor
		
		int i;
		nRows = K.nRows;
		rows = new R2DSparseRow*[nRows];
		indexArray = new int[nRows+1];
		for (i=0; i<nRows;i++) 
		{			
			rows[i] = new R2DSparseRow(*K.rows[i]);
			indexArray[i+1] = -1;
			//*K.indexArray[i+1];	
		}
	}

	 R2DSparseRow::R2DSparseRow (const R2DSparseRow &currentRow) {
	
		int currentRowSize;
		int j;
		int currentColumnIndex;
		double currentValue;

		//currentRowSize=currentRow.entries.size();
		rowNumber=currentRow.rowNumber;
		entries=currentRow.entries;
		currentRowSize=entries.size();
		for (j=0;j<currentRowSize;j++)
			entries[j].P=0;
	}// end of copy constructor

	void R2DSparseMatrix::preconditionedCGs(double *X, double *B)
	{

		double *Q, *D, *R, *Min, *S;
		int i, imax, j;
		double deltaNew = 0, deltaOld = 0, delta0, err, err2delta0, alpha, beta, dTq;

		Q = D = R = S = Min = NULL;

		Q = new double[nRows];
		D = new double[nRows];
		R = new double[nRows];
		Min = new double[nRows];
		S = new double[nRows];

		//initializing variables for the iterative solution
	
		i = 0;
		imax = 2*nRows;	//theortically this method should converge in n iterations or less 
		err = 0.000000001;
		vectorProduct(Q,X);										// {Q} = [A]{X}
		for(j = 0; j < nRows; j++) Min[j] = 1.0/getValue(j+1,j+1); // Preconditioner M[i] = K[i,i]
		for (j = 0; j < nRows; j++) 
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
			vectorProduct(Q,D);									// {Q} = [A]{D}
			dTq = 0;											// sum of {D[i]}*{Q[i]}		
			for(j = 0; j < nRows; j++) dTq = dTq + D[j]*Q[j];
			alpha = deltaNew / dTq;
			for(j = 0; j < nRows; j++) X[j] += alpha*D[j];			// {X} = {X} + alpha*{D}
			if ((i % 50) == 0)									// if divisable by 50
			{
				vectorProduct(Q,X);								// {Q} = [A]{X}
				for(j = 0; j < nRows; j++) R[j] = B[j] - Q[j];		// {R} = {B} - {Q}
			}
			else
				for(j = 0; j < nRows; j++) R[j] -= alpha*Q[j];		// {R} = {R} - alpha*{Q}
			for(j = 0; j < nRows; j++) S[j] = Min[j]*R[j];			// {S} = [Min]{R}
			deltaOld = deltaNew;
			deltaNew = 0;
			for(j = 0; j < nRows; j ++) deltaNew += R[j]*S[j];
			beta = deltaNew / deltaOld;
			for(j = 0; j < nRows; j ++) D[j] = S[j] + beta*D[j];	 // {D} = {S} + beta*{D}
			i++;
		}

	//	CString check;
	//	check.Format("number of iterations = %d", i);
	//	AfxMessageBox(check);

		delete[] Q;
		delete[] D;
		delete[] R;
		delete[] Min;
		delete[] S;

	}

	void R2DSparseMatrix::conjugateGradients(double *X, double *B)
	{
		// the conjugate gradients solution code

		double *Q, *D, *R;
		int i, imax, j;
		double deltaNew = 0, deltaOld = 0, delta0, err, err2delta0, alpha, beta, dTq;

	
		Q = D = R = NULL;

		Q = new double[nRows];
		D = new double[nRows];
		R = new double[nRows];

		//initializing variables for the iterative solution
	
		i = 0;
		imax = 2*nRows;		//theortically this method should converge in n iterations or less 
		err = 0.0000001;
		vectorProduct(Q,X); // {Q} = [A]{X}
		for (j = 0; j < nRows; j++) 
		{
			R[j] = B[j] - Q[j];
			D[j] = R[j];
			deltaNew += R[j]*R[j];
		}
		delta0 = deltaNew;
	//	err2delta0 = err*err*delta0;
		err2delta0 = 0.0000001;

		//the iterative solution

		while ( (i < imax) && (deltaNew > err2delta0) )
		{
			vectorProduct(Q,D);								// {Q} = [A]{D}
			dTq = 0;										// sum of {D[i]}*{Q[i]}		
			for(j = 0; j < nRows; j++) dTq = dTq + D[j]*Q[j];
			alpha = deltaNew / dTq;
			for(j = 0; j < nRows; j++) X[j] += alpha*D[j];		// {X} = {X} + alpha*{D}
			if ((i % 50) == 0)								// if divisable by 50
			{
				vectorProduct(Q,X);						// {Q} = [A]{X}
				for(j = 0; j < nRows; j++) R[j] = B[j] - Q[j];  // {R} = {B} - {Q}
			}
			else
				for(j = 0; j < nRows; j++) R[j] -= alpha*Q[j];  // {R} = {R} - alpha*{Q}
			deltaOld = deltaNew;
			deltaNew = 0;
			for(j = 0; j < nRows; j ++) deltaNew += R[j]*R[j];
			beta = deltaNew / deltaOld;
			for(j = 0; j < nRows; j ++) D[j] = R[j] + beta*D[j]; // {D} = {R} + beta*{D}
			i++;
		}

	//	CString check;
	//	check.Format("number of iterations = %d", i);
	//	AfxMessageBox(check);
	
		delete[] Q;
		delete[] D;
		delete[] R;
	}

	void R2DSparseMatrix::gmres(double *X, double *B, int m, int k, double rec)
	{
		double *Q, *D, *R, *V, *H;
		int i, j, nk, i0, im, it, ii;
		double tem = 1, res = 0, ccos, ssin, res0, resCheck;

		Q = D = R = V = H = NULL;

		Q = new double[nRows];
		D = new double[nRows];
		R = new double[m+1];
		V = new double[nRows*m];
		H = new double[m*(m+1)];

		vectorProduct(Q,X);	 // {Q} = [A]{X}
		
		for(ii=0; ii < nRows; ii++)
		{
			Q[ii] = B[ii] - Q[ii];
			res += Q[ii]*Q[ii];
		}

		res = sqrt(res);
		res0 = res;

		//the iterative solution

		for (it = 0; it < k; it++)
		{
			nk = m;

			//Arnoldi method (which uses the Gram-Schmidt orthogonalization method)
			for(ii=0; ii<nRows; ii++)	Q[ii] /= res;
			for(ii=0; ii<=nk; ii++)		R[ii] = 0;
			R[0] = res;
		
			for(j = 0; j < m; j++)
			{
				for(ii = 0; ii < nRows; ii++) V[ii*m + j] = Q[ii];  
				vectorProduct(D,Q);
				for(ii = 0; ii < nRows; ii++) Q[ii] = D[ii];
				double sum = 0.0; 
				for(ii = 0; ii < nRows; ii++) sum += Q[ii];

				for(i = 0; i <= j; i++)
				{
					H[i*m+j] = 0;
					for(ii = 0; ii < nRows; ii++) H[i*m + j] += Q[ii]*V[ii*m + i];
					for(ii = 0; ii < nRows; ii++) Q[ii] -= H[i*m + j]*V[ii*m + i];
				}
				tem = 0;
				for(ii = 0; ii < nRows; ii++) tem += Q[ii]*Q[ii];
				tem = sqrt(tem);
				H[(j+1)*m + j] = tem;
				if (tem < rec)
				{
					nk = j + 1;
					goto l5;
				}
				for(ii = 0; ii < nRows; ii++) Q[ii] /= tem;
			}
			//triangularization
	l5:		for (i = 0; i < nk; i++)
			{
				im = i + 1;
				tem = (1.0)/(sqrt(H[i*m + i]*H[i*m + i] + H[im*m + i]*H[im*m + i]));
				ccos = H[i*m + i]*tem;
				ssin = -H[im*m + i]*tem;
				for(j = i; j < nk;j++)
				{
					tem = H[i*m + j];
					H[i*m + j] = ccos * tem - ssin * H[im*m + j];
					H[im*m + j] = ssin * tem + ccos *H[im*m + j];
				}
				R[im] = ssin * R[i];
				R[i] *= ccos;
			}
			double sum = 0.0;
			for(ii = 0; ii < (m*(m+1)); ii++) sum += H[ii];
			sum = 0.0; 
			for(ii = 0; ii < (m+1); ii++) sum += R[ii];
			//solution of linear system
			for (i = (nk - 1); i >= 0; i--)
			{
				R[i] /= H[i*m + i];
				for (i0 = (i - 1); i0 >= 0; i0 --) R[i0] -=H[i0*m + i] * R[i];
			}
			sum = 0.0; 
			for(ii = 0; ii < (m+1); ii++) sum += R[ii];
			for (i = 0; i < nk; i++)
			{
				for (ii =0; ii < nRows; ii++) X[ii] += R[i] * V[ii*m + i];	
			}

			//new residual and stopping tests
			vectorProduct(Q,X);	
			res = 0;
			for(ii=0; ii<nRows; ii++)
			{
				Q[ii] = B[ii] - Q[ii];
				res += Q[ii]*Q[ii];
			}
			res = sqrt(res);
			resCheck = res/res0;
		//	if (res < rec) break;
			if(resCheck < rec) break;

		}
	
		
		delete[] Q;
		delete[] D;
		delete[] R;
		delete[] V;
		delete[] H;

	}

	int R2DSparseMatrix::PCJAC_GMRES(double *X, double *B, int m, int k, double rec)
	{

		double *Q, *D, *R, *Min, *V, *H;
		double tem = 1, res = 0, res0, ccos, ssin, resCheck, convergence;
		int ii, it, i, j, nk, i0, im;

		Q = D = R = Min = V = H = NULL;

		Q = new double[nRows];
		D = new double[nRows];
		R = new double[m+1];
		Min = new double[nRows];
		V = new double[nRows*m];
		H = new double[m*(m+1)];

		//initializing variables for the iterative solution

		vectorProduct(Q,X);	 // {Q} = [A]{X}
		for(j = 0; j < nRows; j++) Min[j] = 1.0/getValue(j+1,j+1); // Preconditioner M[i] = K[i,i]

		for (ii = 0; ii < nRows; ii++) 
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
			for (ii = 0; ii < nRows; ii++) Q[ii] = Q[ii] / res;
			for (ii = 0; ii <= nk; ii++) R[ii] = 0;
			R[0] = res;

			for(j = 0; j < m; j ++)
			{
				for (ii = 0; ii < nRows; ii++) V[ii*m + j] = Q[ii];
				vectorProduct(D,Q);		
				for (ii = 0; ii < nRows; ii++) Q[ii] = Min[ii]*D[ii];
		
				for(i = 0; i <= j; i++)
				{
					H[i*m + j] = 0;
					for(ii = 0; ii < nRows; ii++) H[i*m + j] += Q[ii]*V[ii*m + i];
					for(ii = 0; ii < nRows; ii++) Q[ii] -= H[i*m + j]*V[ii*m + i];
				}
				tem = 0;
				for (ii = 0; ii < nRows; ii++) tem += Q[ii]*Q[ii];
				tem = sqrt(tem);
				H[(j+1)*m + j] = tem;
				if (tem < rec)
				{
					nk = j + 1;
					goto l5;
				}
				for (ii = 0; ii < nRows; ii++) Q[ii] = Q[ii] / tem;
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
				for (i0 = i - 1; i0 >= 0; i0--) R[i0] = R[i0] - H[i0*m + i]*R[i];
			}
			for (i = 0; i < nk; i++)
			{
				for(ii = 0; ii < nRows; ii++) 
					X[ii] += R[i]*V[ii*m + i];
			}
			//new residual and stopping tests
			vectorProduct(Q,X);								// {Q} = [A]{X}
			res = 0;
			for(ii = 0; ii < nRows; ii++)
			{
				D[ii] = B[ii] - Q[ii];
				Q[ii] = Min[ii]*D[ii];
				res += Q[ii]*Q[ii];
			}
			res = sqrt(res);
			resCheck = res/res0;

			if(resCheck < rec) break;
			else if (it == (k-1)) ;

		}

		delete[] Q;
		delete[] D;
		delete[] R;
		delete[] Min;
		delete[] V;
		delete[] H;

		return(it+1);

	}

	int R2DSparseMatrix::PCJAC_GMRES(double *B, int m, int k, double rec)
	{

		double *X, *Q, *D, *R, *Min, *V, *H;
		double tem = 1, res = 0, res0, ccos, ssin, resCheck, convergence;
		int ii, it, i, j, nk, i0, im;

		X = Q = D = R = Min = V = H = NULL;

		X = new double[nRows];
		for (ii = 0; ii < nRows; ii++) X[ii] = 0.0;
		Q = new double[nRows];
		D = new double[nRows];
		R = new double[m+1];
		Min = new double[nRows];
		V = new double[nRows*m];
		H = new double[m*(m+1)];

		//initializing variables for the iterative solution

		vectorProduct(Q,X);	 // {Q} = [A]{X}
		for(j = 0; j < nRows; j++) Min[j] = 1.0/getValue(j+1,j+1); // Preconditioner M[i] = K[i,i]

		for (ii = 0; ii < nRows; ii++) 
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
			for (ii = 0; ii < nRows; ii++) Q[ii] = Q[ii] / res;
			for (ii = 0; ii <= nk; ii++) R[ii] = 0;
			R[0] = res;

			for(j = 0; j < m; j ++)
			{
				for (ii = 0; ii < nRows; ii++) V[ii*m + j] = Q[ii];
				vectorProduct(D,Q);		
				for (ii = 0; ii < nRows; ii++) Q[ii] = Min[ii]*D[ii];
		
				for(i = 0; i <= j; i++)
				{
					H[i*m + j] = 0;
					for(ii = 0; ii < nRows; ii++) H[i*m + j] += Q[ii]*V[ii*m + i];
					for(ii = 0; ii < nRows; ii++) Q[ii] -= H[i*m + j]*V[ii*m + i];
				}
				tem = 0;
				for (ii = 0; ii < nRows; ii++) tem += Q[ii]*Q[ii];
				tem = sqrt(tem);
				H[(j+1)*m + j] = tem;
				if (tem < rec)
				{
					nk = j + 1;
					goto l5;
				}
				for (ii = 0; ii < nRows; ii++) Q[ii] = Q[ii] / tem;
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
				for (i0 = i - 1; i0 >= 0; i0--) R[i0] = R[i0] - H[i0*m + i]*R[i];
			}
			for (i = 0; i < nk; i++)
			{
				for(ii = 0; ii < nRows; ii++) 
					X[ii] += R[i]*V[ii*m + i];
			}
			//new residual and stopping tests
			vectorProduct(Q,X);								// {Q} = [A]{X}
			res = 0;
			for(ii = 0; ii < nRows; ii++)
			{
				D[ii] = B[ii] - Q[ii];
				Q[ii] = Min[ii]*D[ii];
				res += Q[ii]*Q[ii];
			}
			res = sqrt(res);
			resCheck = res/res0;

			if(resCheck < rec) break;
			else if (it == (k-1)) ;

		}

		for(ii = 0; ii < nRows; ii++)  B[ii] = X[ii];

		delete[] X;
		delete[] Q;
		delete[] D;
		delete[] R;
		delete[] Min;
		delete[] V;
		delete[] H;

		return(it+1);

	}

	int R2DSparseMatrix::PCILU_GMRES(double *X, double *B, int m, int k, double rec)
	{

		double *Q, *D, *R, *Y, *V, *H;
		double tem = 1, res = 0, res0, ccos, ssin, resCheck, convergence;
		int ii, it, i, j, nk, i0, im;

		Q = D = R = Y = V = H = NULL;

		Q = new double[nRows];
		D = new double[nRows];
		R = new double[m+1];
		Y = new double[nRows];
		V = new double[nRows*m];
		H = new double[m*(m+1)];

		R2DSparseMatrix *ILUmat = NULL;
		ILUmat = new R2DSparseMatrix(*this,nRows);
		ILUmat->LUDecompose();


		//initializing variables for the iterative solution

		vectorProduct(Q,X);	 // {Q} = [A]{X}

		for (ii = 0; ii < nRows; ii++) 
		{
			D[ii] = B[ii] - Q[ii];
			Q[ii] = D[ii];
			res += Q[ii]*Q[ii];
		}

		res = sqrt(res);
		res0 = res;

		//the iterative solution

		for(it = 0; it < k; it++)
		{
			nk = m;
		
			//Arnoldi method (which uses the Gram-Schmidt orthogonalization method)
			for (ii = 0; ii < nRows; ii++) Q[ii] = Q[ii] / res;
			for (ii = 0; ii < nRows; ii++) Y[ii] = 0;
			for (ii = 0; ii <= nk; ii++) R[ii] = 0;
			R[0] = res;

			for(j = 0; j < m; j ++)
			{
				for (ii = 0; ii < nRows; ii++) V[ii*m + j] = Q[ii];
				ILUmat->forwardSweep(Q);
				ILUmat->backSubstitute(Q);
				vectorProduct(D,Q);		
				for (ii = 0; ii < nRows; ii++) Q[ii] = D[ii];
		
				for(i = 0; i <= j; i++)
				{
					H[i*m + j] = 0;
					for(ii = 0; ii < nRows; ii++) H[i*m + j] += Q[ii]*V[ii*m + i];
					for(ii = 0; ii < nRows; ii++) Q[ii] -= H[i*m + j]*V[ii*m + i];
				}
				tem = 0;
				for (ii = 0; ii < nRows; ii++) tem += Q[ii]*Q[ii];
				tem = sqrt(tem);
				H[(j+1)*m + j] = tem;
				if (tem < rec)
				{
					nk = j + 1;
					goto l5;
				}
				for (ii = 0; ii < nRows; ii++) Q[ii] = Q[ii] / tem;
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
				for(ii = 0; ii < nRows; ii++) 
					Y[ii] += R[i]*V[ii*m + i];
			}
			//new residual and stopping tests
			ILUmat->forwardSweep(Y);
			ILUmat->backSubstitute(Y);
			
			for(ii = 0; ii < nRows; ii++) X[ii] +=Y[ii];
			vectorProduct(Q,X);								// {Q} = [A]{X}
			res = 0;
			for(ii = 0; ii < nRows; ii++)
			{
				D[ii] = B[ii] - Q[ii];
				Q[ii] = D[ii];
				res += Q[ii]*Q[ii];
			}
			res = sqrt(res);
			resCheck = res/res0;

			if(resCheck < rec) break;
			else if (it == (k-1)) ;

		}

		delete[] Q;
		delete[] D;
		delete[] R;
		delete[] Y;
		delete[] V;
		delete[] H;
		delete ILUmat;

		return(it+1);

	}

	int R2DSparseMatrix::PCILU_GMRES(double *B, int m, int k, double rec)
	{

		double *X, *Q, *D, *R, *Y, *V, *H;
		double tem = 1, res = 0, res0, ccos, ssin, resCheck, convergence;
		int ii, it, i, j, nk, i0, im;

		X = Q = D = R = Y = V = H = NULL;

		X = new double[nRows];
		for (ii = 0; ii < nRows; ii++) X[ii] = 0.0;
		Q = new double[nRows];
		D = new double[nRows];
		R = new double[m+1];
		Y = new double[nRows];
		V = new double[nRows*m];
		H = new double[m*(m+1)];

		R2DSparseMatrix *ILUmat = NULL;
		ILUmat = new R2DSparseMatrix(*this,nRows);
		ILUmat->LUDecompose();

		//initializing variables for the iterative solution

		vectorProduct(Q,X);	 // {Q} = [A]{X}
		
		for (ii = 0; ii < nRows; ii++) 
		{
			D[ii] = B[ii] - Q[ii];
			Q[ii] = D[ii];
			res += Q[ii]*Q[ii];
		}

		res = sqrt(res);
		res0 = res;

		//the iterative solution

		for(it = 0; it < k; it++)
		{
			nk = m;
		
			//Arnoldi method (which uses the Gram-Schmidt orthogonalization method)
			for (ii = 0; ii < nRows; ii++) Q[ii] = Q[ii] / res;
			for (ii = 0; ii < nRows; ii++) Y[ii] = 0;
			for (ii = 0; ii <= nk; ii++) R[ii] = 0;
			R[0] = res;

			for(j = 0; j < m; j ++)
			{
				for (ii = 0; ii < nRows; ii++) V[ii*m + j] = Q[ii];
				ILUmat->forwardSweep(Q);
				ILUmat->backSubstitute(Q);
				vectorProduct(D,Q);		
				for (ii = 0; ii < nRows; ii++) Q[ii] = D[ii];
		
				for(i = 0; i <= j; i++)
				{
					H[i*m + j] = 0;
					for(ii = 0; ii < nRows; ii++) H[i*m + j] += Q[ii]*V[ii*m + i];
					for(ii = 0; ii < nRows; ii++) Q[ii] -= H[i*m + j]*V[ii*m + i];
				}
				tem = 0;
				for (ii = 0; ii < nRows; ii++) tem += Q[ii]*Q[ii];
				tem = sqrt(tem);
				H[(j+1)*m + j] = tem;
				if (tem < rec)
				{
					nk = j + 1;
					goto l5;
				}
				for (ii = 0; ii < nRows; ii++) Q[ii] = Q[ii] / tem;
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
				for(ii = 0; ii < nRows; ii++) 
					Y[ii] += R[i]*V[ii*m + i];
			}
			//new residual and stopping tests
			ILUmat->forwardSweep(Y);
			ILUmat->backSubstitute(Y);
			
			for(ii = 0; ii < nRows; ii++) X[ii] +=Y[ii];
			vectorProduct(Q,X);								// {Q} = [A]{X}
			res = 0;
			for(ii = 0; ii < nRows; ii++)
			{
				D[ii] = B[ii] - Q[ii];
				Q[ii] = D[ii];
				res += Q[ii]*Q[ii];
			}
			res = sqrt(res);
			resCheck = res/res0;

			if(resCheck < rec) break;
			else if (it == (k-1)) ;

		}
		
		for(ii = 0; ii < nRows; ii++)  B[ii] = X[ii];

		delete[] X;
		delete[] Q;
		delete[] D;
		delete[] R;
		delete[] Y;
		delete[] V;
		delete[] H;
		delete ILUmat;

		return(it+1);

	}

	void R2DSparseMatrix::buildProblem(int itnum, int varnum, int elcode, double *RHS)
	 {
		int i, j, k, ii, jj, vectorRow;
		struct belement         *belp;
		struct element			*elp;
		int						ntf, ns;
		
		int nRows, nCols;
		nRows = 3*N.vars;
		nCols = 3*N.vars;
		int *rowIndices = 0;
		int *colIndices = 0;
		double *matrix = 0;
		rowIndices = new int[nRows];
		colIndices = new int[nCols];
		matrix = new double[nRows*nCols];

		elp = gp.El ;
	        
		for(k=0;k<N.elms;k++)
		{
			if(jacobianType == 0)
				ns = get_PGKeAnalJ(elp,elp,varnum,elcode,&ntf,itnum) ;
			else
 				ns = get_PGKeNumJ(elp,elp,varnum,elcode,&ntf,itnum) ;
			for (i = 0; i < elp->nnds; i++)
			{
				for (j = 0; j < N.vars; j++)
				{
					rowIndices[i*N.vars+j] = elp->nps[i]->i * N.vars + (j+1);
					colIndices[i*N.vars+j] = elp->nps[i]->i * N.vars + (j+1);
				}
			}
			for(i = 0; i < ntf; i ++)  
			{
				for (j = 0; j < ns; j++)
				{	
					for (ii = 0; ii < N.vars; ii++)
					{
						for(jj = 0; jj < N.vars ; jj++)
							matrix[(i*N.vars+ii)*nCols + j*N.vars+jj] = SE(i*N.vars+ii, j*N.vars+jj);		
					}
				}
				vectorRow = elp->nps[i]->i * N.vars;
				for (ii = 0; ii < N.vars; ii++)	RHS[vectorRow + ii] += FE(i*N.vars+ii);
			}
			addSubMatrix(nRows,nCols,rowIndices,colIndices,matrix);
			elp = elp->nextelp ;
		}

	delete [] matrix;
	delete [] rowIndices;
	delete [] colIndices;
	matrix = 0;
	colIndices = 0;
	rowIndices = 0;
	nRows = 2*N.vars;
	nCols = 2*N.vars;
	rowIndices = new int[nRows];
	colIndices = new int[nCols];
	matrix = new double[nRows*nCols];

		belp = gp.B ;

		for(k=0;k<N.belms;k++)
		{ 
			if(jacobianType == 0)
				ns = get_bKeAnalJ(belp,belp,varnum,elcode,&ntf,itnum) ;
			else
				ns = get_bKeNumJ(belp,belp,varnum,elcode,&ntf,itnum) ;
			for (i = 0; i < belp->nnds; i++)
			{
				for (j = 0; j < N.vars; j++)
				{
					rowIndices[i*N.vars+j] = belp->nps[i]->i * N.vars + (j+1);
					colIndices[i*N.vars+j] = belp->nps[i]->i * N.vars + (j+1);
				}
			}
			for(i = 0; i < ntf; i ++)
			{
				for (j = 0; j < ns; j++)
				{
					for (ii = 0; ii < N.vars; ii++)
					{
						for(jj = 0; jj < N.vars ; jj++)
							matrix[(i*N.vars+ii)*nCols + j*N.vars+jj] = SE(i*N.vars+ii, j*N.vars+jj);	
					}
				}
				vectorRow = belp->nps[i]->i * N.vars;
				for (ii = 0; ii < N.vars; ii++)	RHS[vectorRow + ii] += FE(i*N.vars+ii);
			}
			addSubMatrix(nRows,nCols,rowIndices,colIndices,matrix);
			belp = belp->nextbelp ;
		}

	delete [] matrix;
	delete [] rowIndices;
	delete [] colIndices;
	matrix = 0;
	colIndices = 0;
	rowIndices = 0;

	// code for looking at stats for sparse matrix

	/*	int nRows = getNRows();
		R2DSparseRow *currentRow;
		int currentRowSize;
		int maxRowSize = 0, averageRowSize = 0, sumOfRows = 0;
		int minRowSize;
		currentRow = getRow(1);
		minRowSize = currentRow->getRowSize();
		int count[50];

		for(i=0;i<50;i++)count[i]=0;
		
		for(i=0;i<nRows;i++)
		{
			currentRow = getRow(i+1);
			currentRowSize = currentRow->getRowSize();
			if (currentRowSize > maxRowSize)
				maxRowSize = currentRowSize;
			if (currentRowSize < minRowSize)
				minRowSize = currentRowSize;
			sumOfRows += currentRowSize;
			count[currentRowSize] +=1;

		}
		averageRowSize = sumOfRows / nRows;

			CString check;
			check.Format("max %d\t  min %d\t average %d\t numrows %d\t",maxRowSize,minRowSize,averageRowSize, nRows);
			AfxMessageBox(check);

			ofstream textfile ("histogram.txt");
			for (i=0;i<50;i++) textfile << count[i] << "\n";
			textfile.close();*/
	 }

	int PCILU_gmres(int nset, int itnum, int varnum, int elcode, int m, int k, double rec/*, FILE* fp*/)
	{
		int i, j, nk, i0, im, it, ii;
		double tem = 1, res=0, res0, resCheck, zero, ccos, ssin, eps, convergence;
		double	varsum,oldH;
		double	*dp;

		R2DSparseMatrix *KPrime = NULL;
		R2DSparseMatrix *ILUmat = NULL;
	
		KPrime = new R2DSparseMatrix(N.vars*N.nodes, 21);
	
		zero = rec;
		eps = 1.e-5;    
					
		KPrime->buildProblem(itnum,varnum,1,eqnsets[nset].Bp);
		
		ILUmat=new R2DSparseMatrix(*KPrime,N.vars*N.nodes); // invokes copy constructor
			
		ILUmat->LUDecompose();

		for(ii=0; ii<Nukns; ii++)	eqnsets[nset].Fp[ii] = 0.0;
			
		KPrime->matrixVectorMultSparse(eqnsets[nset].Fp,eqnsets[nset].Qp);

		for(ii=0; ii<Nukns; ii++)
		{
		
			eqnsets[nset].Dp[ii] = eqnsets[nset].Bp[ii] - eqnsets[nset].Qp[ii];//residual vector
		
			eqnsets[nset].Qp[ii] = eqnsets[nset].Dp[ii];
			res += eqnsets[nset].Qp[ii]*eqnsets[nset].Qp[ii];
		}
		res = sqrt(res);
		res0 = res;

		//loop that governs number of reinitializations k
		//could eventually be a while statement...

	

		for (it = 0; it < k; it++)
		{
			nk = m;

			//Arnoldi method (which uses the Gram-Schmidt orthogonalization method)
			for(ii=0; ii<Nukns; ii++)	eqnsets[nset].Qp[ii] /= res;
			for(ii=0; ii<Nukns; ii++)	eqnsets[nset].Yp[ii] = 0;
			for(ii=0; ii<=nk; ii++)		eqnsets[nset].Rp[ii] = 0;
			eqnsets[nset].Rp[0] = res;
		
			for(j = 0; j < m; j++)
			{
				for(ii = 0; ii < Nukns; ii++) Vp(ii, j) = eqnsets[nset].Qp[ii];  
			
				ILUmat->forwardSweep (eqnsets[nset].Qp);
				ILUmat->backSubstitute (eqnsets[nset].Qp);
	
				KPrime->matrixVectorMultSparse(eqnsets[nset].Qp,eqnsets[nset].Dp);
			
				for(ii = 0; ii < Nukns; ii++) eqnsets[nset].Qp[ii] = eqnsets[nset].Dp[ii];
				for(i = 0; i <= j; i++)
				{
					Hp(i,j) = 0;
					for(ii = 0; ii < Nukns; ii++)	Hp(i,j) += eqnsets[nset].Qp[ii] * Vp(ii, i);
					for(ii = 0; ii < Nukns; ii++)	eqnsets[nset].Qp[ii] -= Hp(i,j) * Vp(ii, i);
				}
				tem = 0;
				for(ii = 0; ii < Nukns; ii++) tem += eqnsets[nset].Qp[ii]*eqnsets[nset].Qp[ii];
				tem = sqrt(tem);
				Hp(j+1, j) = tem;
				if (tem < rec)
				{
					nk = j + 1;
					goto l5;
				}
				for(ii = 0; ii < Nukns; ii++) eqnsets[nset].Qp[ii] /= tem;
			}
			//triangularization
	l5:		for (i = 0; i < nk; i++)
			{
				im = i + 1;
				tem = (1.0)/( sqrt(Hp(i,i)*Hp(i,i) + Hp(im,i)*Hp(im,i)) );
				ccos = Hp(i,i)*tem;
				ssin = -Hp(im,i)*tem;
				for(j = i; j < nk;j++)
				{
					tem = Hp(i,j);
					Hp(i,j) = ccos * tem - ssin * Hp(im,j);
					Hp(im,j) = ssin * tem + ccos * Hp(im,j);
				}
				eqnsets[nset].Rp[im] = ssin * eqnsets[nset].Rp[i];
				eqnsets[nset].Rp[i] *= ccos;
			}
			//solution of linear system
			for (i = (nk - 1); i >= 0; i--)
			{
				eqnsets[nset].Rp[i] /= Hp(i,i);
				for (i0 = i -1; i0 >= 0; i0 --)
					eqnsets[nset].Rp[i0] -=Hp(i0,i) * eqnsets[nset].Rp[i];
			}
			//for (ii =0; ii < Nukns; ii++) eqnsets[nset].Qp[ii]=0;
			for (i = 0; i < nk; i++)
			{
			
				for (ii =0; ii < Nukns; ii++) eqnsets[nset].Yp[ii] += eqnsets[nset].Rp[i] * Vp(ii,i);	
			}
			ILUmat->forwardSweep (eqnsets[nset].Yp);
			ILUmat->backSubstitute (eqnsets[nset].Yp);

			for (ii =0; ii < Nukns; ii++) eqnsets[nset].Fp[ii] += eqnsets[nset].Yp[ii];

			//new residual and stopping tests
			KPrime->matrixVectorMultSparse(eqnsets[nset].Fp,eqnsets[nset].Qp);
			//matrixVectorMult(eqnsets[nset].Fp, eqnsets[nset].Qp, itnum, varnum, -1);
			res = 0;
			for(ii=0; ii<Nukns; ii++)
			{
				eqnsets[nset].Dp[ii] = eqnsets[nset].Bp[ii] - eqnsets[nset].Qp[ii];
				eqnsets[nset].Qp[ii] = eqnsets[nset].Dp[ii];
				res += eqnsets[nset].Qp[ii]*eqnsets[nset].Qp[ii];
			}
			res = sqrt(res);
			resCheck = res/res0;
			convergence = 1.0 - (res/res0);
		
		
		
	//		fprintf(fp, " %d\t %f\t %f\t %f\t %f\t %f\n", it+1, res, res0, resCheck, convergence, rec); 
			if (resCheck < rec) break;
			else if (it == (k-1)) ;//	printf("Stopped before convergence\n");
		}
	
		delete KPrime;
		delete ILUmat;

	

	//	CString check;
	//	check.Format("%f",duration);
	//	AfxMessageBox(check);

	//	fprintf(fp," %f\t %f\t %f\t %f\t %d\t %d \t %f\t", rec, res0, res, convergence, it, nk, duration);
	
		// update nodal values based on delta phi (Fp) and set uchange

		dp = eqnsets[nset].Fp ;
	
		for(i=0;i<N.nodes;i++) 
		{
			for(j=0;j<N.vars;j++) 
			{
				gp.iptrs[i]->u[j] += *dp ;
				dp++;
			}	
		}
		
		return(it+1);
	}

	int PCJAC_gmres(int nset, int itnum, int varnum, int elcode, int m, int k, double rec/*, FILE* fp*/)
	{
		int i, j, nk, i0, im, it, ii;
		double tem = 1, res=0, res0, resCheck, zero, ccos, ssin, eps, convergence;
		double	varsum;
		double	*dp;
		//clock_t start, finish;
		double duration;

		R2DSparseMatrix *KPrime = NULL;
		KPrime = new R2DSparseMatrix(N.vars*N.nodes, 21);
	
		zero = rec;
		eps = 1.e-5;

		KPrime->buildProblem(itnum,varnum,1,eqnsets[nset].Bp);
		KPrime->matrixVectorMultSparse(eqnsets[nset].Fp,eqnsets[nset].Qp);
		JacobiPC(eqnsets[nset].Mp, itnum, varnum, -1);

		for(ii=0; ii<Nukns; ii++)
		{
			eqnsets[nset].Dp[ii] = eqnsets[nset].Bp[ii] - eqnsets[nset].Qp[ii];
			eqnsets[nset].Qp[ii] = eqnsets[nset].Mp[ii]*eqnsets[nset].Dp[ii];
			res += eqnsets[nset].Qp[ii]*eqnsets[nset].Qp[ii];
		}
		res = sqrt(res);
		res0 = res;

		//loop that governs number of reinitializations k
		//could eventually be a while statement...

		//start = clock();

		for (it = 0; it < k; it++)
		{
			nk = m;

			//Arnoldi method (which uses the Gram-Schmidt orthogonalization method)
			for(ii=0; ii<Nukns; ii++)	eqnsets[nset].Qp[ii] /= res;
			for(ii=0; ii<=nk; ii++)		eqnsets[nset].Rp[ii] = 0;
			eqnsets[nset].Rp[0] = res;
		
			for(j = 0; j < m; j++)
			{
				for(ii = 0; ii < Nukns; ii++) Vp(ii, j) = eqnsets[nset].Qp[ii];  
				KPrime->matrixVectorMultSparse(eqnsets[nset].Qp,eqnsets[nset].Dp);
				//matrixVectorMult(eqnsets[nset].Qp,eqnsets[nset].Dp, itnum, varnum, -1);
				for(ii = 0; ii < Nukns; ii++) eqnsets[nset].Qp[ii] = eqnsets[nset].Mp[ii]*eqnsets[nset].Dp[ii];

				for(i = 0; i <= j; i++)
				{
					Hp(i,j) = 0;
					for(ii = 0; ii < Nukns; ii++)	Hp(i,j) += eqnsets[nset].Qp[ii] * Vp(ii, i);
					for(ii = 0; ii < Nukns; ii++)	eqnsets[nset].Qp[ii] -= Hp(i,j) * Vp(ii, i);
				}
				tem = 0;
				for(ii = 0; ii < Nukns; ii++) tem += eqnsets[nset].Qp[ii]*eqnsets[nset].Qp[ii];
				tem = sqrt(tem);
				Hp(j+1, j) = tem;
				if (tem < rec)
				{
					nk = j + 1;
					goto l5;
				}
				for(ii = 0; ii < Nukns; ii++) eqnsets[nset].Qp[ii] /= tem;
			}
			//triangularization
	l5:		for (i = 0; i < nk; i++)
			{
				im = i + 1;
				tem = (1.0)/( sqrt(Hp(i,i)*Hp(i,i) + Hp(im,i)*Hp(im,i)) );
				ccos = Hp(i,i)*tem;
				ssin = -Hp(im,i)*tem;
				for(j = i; j < nk;j++)
				{
					tem = Hp(i,j);
					Hp(i,j) = ccos * tem - ssin * Hp(im,j);
					Hp(im,j) = ssin * tem + ccos * Hp(im,j);
				}
				eqnsets[nset].Rp[im] = ssin * eqnsets[nset].Rp[i];
				eqnsets[nset].Rp[i] *= ccos;
			}
			//solution of linear system
			for (i = (nk - 1); i >= 0; i--)
			{
				eqnsets[nset].Rp[i] /= Hp(i,i);
				for (i0 = i -1; i0 >= 0; i0 --)
					eqnsets[nset].Rp[i0] -=Hp(i0,i) * eqnsets[nset].Rp[i];
			}
			for (i = 0; i < nk; i++)
			{
				for (ii =0; ii < Nukns; ii++) eqnsets[nset].Fp[ii] += eqnsets[nset].Rp[i] * Vp(ii,i);	
			}

			//new residual and stopping tests
			KPrime->matrixVectorMultSparse(eqnsets[nset].Fp,eqnsets[nset].Qp);
			res = 0;
			for(ii=0; ii<Nukns; ii++)
			{
				eqnsets[nset].Dp[ii] = eqnsets[nset].Bp[ii] - eqnsets[nset].Qp[ii];
				eqnsets[nset].Qp[ii] = eqnsets[nset].Mp[ii]*eqnsets[nset].Dp[ii];
				res += eqnsets[nset].Qp[ii]*eqnsets[nset].Qp[ii];
			}
			res = sqrt(res);
			resCheck = res/res0;
			convergence = 1.0 - (res/res0);
	//		fprintf(fp, " %d\t %f\t %f\t %f\t %f\t %f\n", it+1, res, res0, resCheck, convergence, rec); 
			if (resCheck < rec) break;
			else if (it == (k-1)) ;//	printf("Stopped before convergence\n");
		}
		delete KPrime;
	
		//finish = clock();
		//duration = (double)(finish - start) / CLOCKS_PER_SEC;

	//	CString check;
	//	check.Format("%f",duration);
	//	AfxMessageBox(check);

	//	fprintf(fp," %f\t %f\t %f\t %f\t %d\t %d \t %f\t", rec, res0, res, convergence, it, nk, duration);
	
		// update nodal values based on delta phi (Fp) and set uchange

		dp = eqnsets[nset].Fp ;
		for(i=0;i<N.nodes;i++){
			for(j=0;j<N.vars;j++){
				gp.iptrs[i]->u[j] += *dp ;
				dp++;
			}
		}
	


		return(it+1);
	}
}

