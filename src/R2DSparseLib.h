//
//		R2DSparseLib.h
//
//		created:	Sept. 23, 2006
//		revised:	Sept. 23, 2006
//
//		This library implements sparse matrix methods for solution of River2D problems.
//

#ifndef R2DSPARSELIB
#define R2DSPARSELIB

#include <vector>
#include "HabitatTIN.h"
//#include "ElementMtx.h"

namespace R2DSparseLib {

	struct R2DSparseValue {
		double value;
		int columnIndex;
		int P;
	};

	class R2DSparseRow {
	private:
		int rowNumber;
		std::vector<R2DSparseValue> entries;
	public:
		R2DSparseRow(int row, int initialSize, int maxColNumber);
		~R2DSparseRow();
		int getRowNumber() {return rowNumber;};
		void setRowNumber(int n) {rowNumber = n;};
		int getRowSize() {return (int) entries.size();};
		int getRowCapacity() {return (int) entries.capacity();};
		void setEntry (int index, int ncolumn, double value);
		void setZeroEntry(int ncolumn);
		void setValue(int ncolumn, double value);
		double getValue(int ncolumn);
		void addValue(int ncolumn, double value);
		void addRow(double factor, R2DSparseRow *pivotRow, int *indexArray );
		void reduceRow(R2DSparseRow *pivotRow, int *indexArray );
		void sweepRow(double *rhs);
		void backSubRow(double *rhs);
		double vectorProduct(double *vector);
		double getValueByIndex(int index){return entries[index].value;};
		int getColumnByIndex(int index) {return entries[index].columnIndex;};
		int getPbyIndex(int index) {return entries[index].P;};
		void clearRow();
		double getValueByColumn(int column);
		R2DSparseRow::R2DSparseRow (const R2DSparseRow &currentRow); 
		void testRow();
		
	};

	class R2DSparseMatrix {
	private:
		int	nRows;
		R2DSparseRow ** rows;
		int *indexArray;
	public:
		R2DSparseMatrix(int numRows, int rowSize);
		~R2DSparseMatrix();
		int	getNRows() {return nRows;};
		int	getNEntries();
		int *getIndexArray() {return indexArray;}
		R2DSparseRow* getRow(int i) {return rows[i-1];}
		void setValue(int row, int column, double value);
		double getValue(int row, int column);
		void addValue(int row, int column, double value);
		void addRow(int pivotRow, int targetRow, double factor);
		void addSubMatrix(int nr, int nc, const int *rownums, const int *colnums, double *sMatrix);
		//void addSubMatrix(ElementMatrix* sMatrix);
		void reduceRow(int pivotRow, int targetRow);
		void LUDecompose();
		void forwardSweep(double *rhs);
		void backSubstitute(double *rhs);
		void vectorProduct(double *result, double *vector);
		void buildProblem(int itnum, int varnum, int elcode, double *RHS);
		void matrixVectorMultSparse(double *X, double *Q);
		void preconditionedCGs(double *X, double *B);
		void conjugateGradients(double *X, double *B);
		int PCJAC_GMRES(double *X, double *B, int m, int k, double rec);
		int PCJAC_GMRES(double *B, int m, int k, double rec);
		int PCILU_GMRES(double *X, double *B, int m, int k, double rec);
		int PCILU_GMRES(double *B, int m, int k, double rec);
		void gmres(double *X, double *B, int m, int k, double rec);
		void clearMatrix();

		R2DSparseMatrix::R2DSparseMatrix (const R2DSparseMatrix &K,int numRows);
		
	};

	//compiler directives so that functions call be declared in a .cpp file but called in a .c file
	#ifdef __cplusplus
	 extern "C" {
	 #endif
	 
	 #if defined(__STDC__) || defined(__cplusplus)
	   extern int PCILU_gmres(int nset, int itnum, int varnum, int elcode, int m, int k, double rec/*, FILE* fp*/);
	   extern int PCJAC_gmres(int nset, int itnum, int varnum, int elcode, int m, int k, double rec/*, FILE* fp*/);
	  #else
		 extern int PCILU_gmres();
		 extern int PCJAC_gmres();
	#endif
	 
	 #ifdef __cplusplus
	 }
	 #endif

}

#endif // R2DSPARSELIB
