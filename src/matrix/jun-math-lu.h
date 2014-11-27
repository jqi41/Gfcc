#ifndef _JUN_MATH_LU_H
#define _JUN_MATH_LU_H	1

#include "../matrix/jun-math-lu.h"

namespace jun {

/* LU decomposition */
template<typename T>
bool LUDecompose(Matrix<T> *&A, Matrix<T> *&L, Matrix<T> *&U) {
	if (A->GetColNum() != A->GetRowNum()) {
		cout << "Matrix A should be a squared matrix" << endl;
		return false;
	}
	if (L==0) {
		L = new Matrix<T>;
		L->CreateEyeMat(A->GetRowNum());
	}
	if (U==0) {
		U = new Matrix<T>;
		U->CreateMat(A->GetRowNum(), A->GetColNum());
	}
	for (int k=1; k<=A->GetRowNum()-1; k++) {
		if (A->data[k-1][k-1] == 0) {
			cout << "A is a singular matrix. Have to stop at this position." << endl;
			return false;
		}
		for (int i=k+1; i<=A->GetColNum(); i++) {
			for (int j=1; j<=k-1; j++) {
				A->data[i-1][k-1] = A->data[i-1][k-1] - A->data[j-1][k-1]*A->data[i-1][j-1];
			}
			A->data[i-1][k-1] = A->data[i-1][k-1] / A->data[k-1][k-1];
		}
		for (int j=k+1; j<=A->GetColNum(); j++) {
			for (int i=1; i<=k; i++) {
				A->data[k][j-1] = A->data[k][j-1] - A->data[k][i-1]*A->data[i-1][j-1];
			}
		}
	}
	for (int i=2; i<=A->GetRowNum(); i++) {
		for (int j=1; j<=i-1; j++) {
			L->data[i-1][j-1] = A->data[i-1][j-1];
		}
	}
	for (int i=0; i<A->GetRowNum(); i++) {
		for (int j=i; j<A->GetRowNum(); j++) {
			U->data[i][j] = A->data[i][j];
		}
	}
	return true;
}

/* Forward Iteration */
template<typename T>
bool ForwardIter(Matrix<T> *L, Vector<T> *b, Vector<T> *&y) {
	if (L==0 || b==0) {
		cout << "Matrix L or b should be initialized at first." << endl;
		return false;
	}
	if (y==0)
		y = new Vector<T>;
	y->New(b->size());

	for (int i=1; i<=L->GetRowNum(); i++) {
		y->data[i-1] = b->data[i-1];
		for (int j=1; j<=i-1; j++) {
			y->data[i-1] = y->data[i-1] - L->data[i-1][j-1]*y->data[j-1];
		}
	}
	return true;
}

/* Backward Iteration */
template<typename T>
bool BackIter(Matrix<T> *U, Vector<T> *b, Vector<T> *&x) {
	if (U==0 || b==0) {
		cout << "Matrix U and vector b should be initialized at first." << endl;
		return false;
	}
	if (x==0)
		x = new Vector<T>;
	x->New(b->size());

	int32 i = U->GetRowNum();
	int32 j = i;
	while (i>0) {
		if (U->data[i-1][i-1]==0) {
			cout << "The matrix is a singular matrix, stop at this place." << endl;
			return false;
		}
		x->data[i-1] = b->data[i-1];
		while (j>i) {
			x->data[i-1] = x->data[i-1] - U->data[i-1][j-1] * x->data[j-1];
			j--;
		}
		x->data[i-1] = x->data[i-1] / U->data[i-1][i-1];
		i--;
	}
	return true;
}

/* Solving a linear equation Ax=b */
template<typename T>
void LinSolve(Matrix<T> *&A, Vector<T> *b, Vector<T> *&x) {
	Matrix<float32> *L = 0;
	Matrix<float32> *U = 0;
	Vector<float32> *y = 0;
	LUDecompose<float32>(A, L, U);
	ForwardIter<float32>(L, b, y);
	BackIter<float32>(U, y, x);
	L->deleteMat(); delete L; L=0;
	U->deleteMat(); delete U; U=0;
	y->delVec();	delete y; y=0;
	return;
}

/* Inverting a matrix */
template<typename T>
bool InvMat(Matrix<T> *&mat, Matrix<T> *&inv_mat) {
	if (mat==0) {
		cout << "The input matrix should be initialized at first." << endl;
		return false;
	}
	if (inv_mat==0)
		inv_mat = new Matrix<T>;
	inv_mat->CreateEyeMat(mat->GetRowNum());
	int32 n = mat->GetRowNum();

	for (int k=1; k<=mat->GetRowNum(); k++) {
		if (mat->data[k-1][k-1]==0) {
			cout << "The matrix should be not singular. Stop there." << endl;
			return false;
		}
		for (int j=k+1; j<=n; j++) {
			mat->data[k-1][j-1] = mat->data[k-1][j-1] / mat->data[k-1][k-1];
		}
		for (int j=1; j<=k; j++) {
			inv_mat->data[k-1][j-1] = inv_mat->data[k-1][j-1] / mat->data[k-1][k-1];
		}
		for (int i=1; i<=n; i++) {
			if (i!=k) {
				for (int j=k+1; j<=n; j++) {
					mat->data[i-1][j-1] = mat->data[i-1][j-1] - mat->data[i-1][k-1]*mat->data[k-1][j-1];
				}
				for (int j=1; j<=k; j++) {
					inv_mat->data[i-1][j-1] = inv_mat->data[i-1][j-1] - mat->data[i-1][k-1]*inv_mat->data[k-1][j-1];
				}
			}
		}
	}
	return true;
}
}

#endif 