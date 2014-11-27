#ifndef _JUN_MATH_H
#define _JUN_MATH_H	1

#include "../matrix/jun-math.h"

namespace jun {

/* Givens rotation (left) */
template<typename T>
void mxGivens_col(Vector<T> *&vec, Matrix<T> *&H, int32 dim_ref, int32 dim_eli) {
	int32 dim = vec->size();
	if (H==0) 
		H = new Matrix<T>;
	H->CreateMat(dim, dim);
	for (int i=0; i<dim; i++) 
		H->data[i][i] = 1.0;
	float32	norm = (vec->data[dim_ref]*vec->data[dim_ref] + vec->data[dim_eli]*vec->data[dim_eli]);
	norm = sqrt(norm);
	float32 c = vec->data[dim_ref] / norm;
	float32 s = vec->data[dim_eli] / norm;

	H->data[dim_ref][dim_ref] = c;
	H->data[dim_eli][dim_eli] = c;
	H->data[dim_ref][dim_eli] = s;
	H->data[dim_eli][dim_ref] = -s;
	return;
}

/* Givens rotation (right) */
template<typename T>
void mxGivens_row(Vector<T> *&vec, Matrix<T> *&H, int32 dim_ref, int32 dim_eli) {
	Matrix<T> *H_tmp = 0;

	mxGivens_col<T>(vec, H_tmp, dim_ref, dim_eli);
	mxTrMat(H_tmp, H);
	
	H_tmp->deleteMat();	
	delete H_tmp;	H_tmp = 0;
}

/* Householder transformation */
template<typename T>
void mxHholder(Vector<T> *vec, Matrix<T> *&Hholder, int m=0) {
	Vector<T> *di = new Vector<T>;
	di->New(vec->size());	
	if (vec->data[0] > 0)	
		di->data[0] = vec->norm2();
	else
		di->data[0] = -1 * vec->norm2();

	Vector<T> *weights = 0;
	Vector<T> *weights_tmp = new Vector<T>;
	T InMat_tmp = 0.0;
	Matrix<T> *OutMat_tmp = 0;
	
	VectorAdd<T>(vec, di, weights);
	
	weights_tmp->New(m+vec->size());
	for (int i=0; i<vec->size(); i++) {
		weights_tmp->data[i+m] = weights->data[i];
	}

	Matrix<T> *diagMat = new Matrix<T>;
	diagMat->CreateMat(vec->size()+m, vec->size()+m);
	for (int n=0; n<vec->size()+m; n++) 
		diagMat->data[n][n] = 1.0;
	
	mxOutProd<T>(weights_tmp, weights_tmp, OutMat_tmp);
	mxInProd<T>(weights_tmp, weights_tmp, InMat_tmp);
	ConstMulMat<T>(2/InMat_tmp, OutMat_tmp);
	mxMatSubMat<T>(diagMat, OutMat_tmp, Hholder);

	/* delete */
	di->delVec();	delete di;	di = 0;
	OutMat_tmp->deleteMat();	delete OutMat_tmp;	OutMat_tmp = 0;
	weights_tmp->delVec();		delete weights_tmp;	weights_tmp = 0;
	weights->delVec();			delete weights;		weights = 0;
	diagMat->deleteMat();		delete diagMat;		diagMat = 0;

	return;
}

/* Singular value decomposition (stage 1) */
template<typename T>
void svd1(Matrix<T> *&A, Matrix<T> *&biDiag_mat, Matrix<T> *&Hholder_col, Matrix<T> *&Hholder_row) {
	Matrix<T> *Hholder_tmp = 0;
	if (biDiag_mat==0)
		biDiag_mat = new Matrix<T>;
	int32 dim_col = A->GetColNum();
	int32 dim_row = A->GetRowNum();
	biDiag_mat->CreateMat(dim_col, dim_col);

	for (int i=0; i<dim_col; i++) {
		/* Dealing with column vector of A */
		Vector<T> *vec_tmp = new Vector<T>;
		vec_tmp->New(A->GetRowNum()-i);
		for (int n=0; n<vec_tmp->size(); n++)
			vec_tmp->data[n] = A->data[n+i][i];
		mxHholder(vec_tmp, Hholder_tmp, i);
		mxLItrMatMul(A, Hholder_tmp);
		if (Hholder_col==0) {
			Hholder_col = new Matrix<T>;
			Hholder_col->CreateMat(dim_row, dim_row);
			for (int m=0; m<dim_row; m++)
				Hholder_col->data[m][m] = 1.0;
		}
		mxLItrMatMul(Hholder_col, Hholder_tmp);
		Hholder_tmp->deleteMat();
		vec_tmp->delVec();

		/* Dealing with row vector of A */
		int err = dim_col - i;
		if (err > 1) {
			vec_tmp->New(A->GetColNum()-i-1);
			for (int n=i+1; n<A->GetColNum(); n++) 
				vec_tmp->data[n-i-1] = A->data[i][n];
			mxHholder(vec_tmp, Hholder_tmp, i+1);
			mxRItrMatMul(A, Hholder_tmp);
			if (Hholder_row==0) {
				Hholder_row = new Matrix<T>;
				Hholder_row->CreateMat(dim_col, dim_col);
				for (int m=0; m<dim_col; m++)
					Hholder_row->data[m][m] = 1.0;
			}
			mxRItrMatMul(Hholder_row, Hholder_tmp);
			Hholder_tmp->deleteMat();
			vec_tmp->delVec();
		}
	}
	for (int n=0; n<dim_col; n++) 
		for (int m=0; m<dim_col; m++)
			biDiag_mat->data[n][m] = A->data[n][m];
	return;
}

/* Singular value decomposition (stage 2) */
template<typename T>
bool svd2(Matrix<T> *&Diag_mat, Matrix<T> *&Givens_left, Matrix<T> *&Givens_right) {
	if (Diag_mat==0) {
		cout << "The diagonal matrix should be initialized at first." << endl;
		return false;
	}
	float32 threshold = 1e-12;
	int32 dim_row = Diag_mat->GetRowNum();
	int32 raw_row = dim_row;
	int32 dim_col = Diag_mat->GetColNum();
	int32 raw_col = dim_col;
	int32 m = 0;
	Matrix<T> *H_row = 0;
	Matrix<T> *H_col = 0;
	Matrix<T> *H_tmp = new Matrix<T>;
	Matrix<T> *H_DiagMat = new Matrix<T>;
	Matrix<T> *Res_DiagMat = new Matrix<T>;
	Vector<T> *vec_tmp = new Vector<T>;

	H_tmp->CreateMat(dim_row, dim_col);
	Res_DiagMat->CreateMat(raw_row, raw_col);

	if (Givens_left==0)
		Givens_left = new Matrix<T>;
	Givens_left->CreateEyeMat(Diag_mat->GetColNum());

	if (Givens_right==0)
		Givens_right = new Matrix<T>;
	Givens_right->CreateEyeMat(Diag_mat->GetRowNum());

	while (dim_row>1) {
		if (dim_col>0) {
			for (m=0; m<dim_col-1; m++) {
				vec_tmp->New(dim_row);
				for (int i=0; i<dim_col; i++)
					vec_tmp->data[i] = Diag_mat->data[m][i];
				mxGivens_row(vec_tmp, H_row, m, m+1);
				mxRItrMatMul(Diag_mat, H_row);
				if (dim_row < raw_row) {
					CopySubMat(H_tmp, H_row, 1);
					mxRItrMatMul(Givens_right, H_tmp);
					H_tmp->zero();
				}
				else 
					mxRItrMatMul(Givens_right, H_row);
				vec_tmp->zero();
			
				for (int n=0; n<dim_row; n++)
					vec_tmp->data[n] = Diag_mat->data[n][m];
				mxGivens_col(vec_tmp, H_col, m, m+1);
				mxLItrMatMul(Diag_mat, H_col);
				if (dim_row<raw_row) {
					CopySubMat(H_tmp, H_col, 1);
					mxLItrMatMul(Givens_left, H_tmp);
				}
				else 
					mxLItrMatMul(Givens_left, H_col);
				H_tmp->zero();
			}
		}
		m--;
		float32 abs_threshold = Diag_mat->data[m][m+1];
		if (fabs(abs_threshold) < threshold) {
			Res_DiagMat->data[m+1][m+1] = Diag_mat->data[m+1][m+1];
			H_DiagMat->CreateMat(m+1, m+1);
		
			CopySubMat(H_DiagMat, Diag_mat, 0);
			Diag_mat->deleteMat(); Diag_mat=0;
			CopyMatrix<T>(Diag_mat, H_DiagMat);
			dim_row = Diag_mat->GetRowNum();
			dim_col = Diag_mat->GetColNum();
			H_row->deleteMat(); H_row = 0;
			H_col->deleteMat(); H_col = 0;
		}
	}
	Res_DiagMat->data[0][0] = Diag_mat->data[0][0];
	Diag_mat->deleteMat();	Diag_mat=0;
	CopyMatrix<T>(Diag_mat, Res_DiagMat);

	vec_tmp->delVec(); vec_tmp = 0;
	H_DiagMat->deleteMat();	H_DiagMat = 0;
	Res_DiagMat->deleteMat();	Res_DiagMat = 0;
	H_tmp->deleteMat(); H_tmp = 0;
	return true;
}

/* Singular value decomposition (main access) */
template<typename T>
bool svd(Matrix<T> *A, Matrix<T> *&U, Matrix<T> *&D, Matrix<T> *&V) {
	if (A==0) {
		cout << "The maxtrix A should be initialized at first." << endl;
		return false;
	}
	if (D==0) 
		D = new Matrix<T>;
	int32 dim_row = A->GetRowNum();
	int32 dim_col = A->GetColNum();
	bool tp_flag = false;

	Matrix<T> *Hholder_col = 0;
	Matrix<T> *Hholder_row = 0;
	Matrix<T> *Givens_left_tmp = 0;
	Matrix<T> *Givens_left = new Matrix<T>;
	Matrix<T> *Givens_right = 0;
	Matrix<T> *A_tr = 0;
	Matrix<T> *D_tmp = 0;
	Matrix<T> *U_tmp = 0;
	Matrix<T> *V_tmp = 0;

	if (dim_row < dim_col) {
		mxTrMat(A, A_tr);
		tp_flag = true;	
	}

	if (tp_flag==0) {
		svd1<float32>(A, D_tmp, Hholder_col, Hholder_row);
		Givens_left->CreateMat(A->GetRowNum(), A->GetRowNum());
	}
	else {
		svd1<float32>(A_tr, D_tmp, Hholder_col, Hholder_row);
		Givens_left->CreateMat(A_tr->GetRowNum(), A_tr->GetRowNum());
	}
	D->CreateMat(A->GetRowNum(), A->GetColNum());
	svd2<float32>(D_tmp, Givens_left_tmp, Givens_right);
	CopySubMat<float32>(D, D_tmp, 0);
	CopySubMat<float32>(Givens_left, Givens_left_tmp, 1);

	if (tp_flag==1) {
		mxMatMulMat(Givens_left, Hholder_col, V_tmp);
		mxTrMat(V_tmp, V);
		mxMatMulMat(Hholder_row, Givens_right, U);
		V_tmp->deleteMat();	V_tmp = 0;
		A_tr->deleteMat();	A_tr = 0;
	}
	else {
		mxMatMulMat(Givens_left, Hholder_col, U_tmp);
		mxTrMat(U_tmp, U);
		mxMatMulMat(Hholder_row, Givens_right, V);
		U_tmp->deleteMat();	U_tmp = 0;
	}
	Hholder_col->deleteMat();		Hholder_col = 0;
	Hholder_row->deleteMat();		Hholder_row = 0;
	Givens_left->deleteMat();		Givens_left = 0;
	Givens_right->deleteMat();		Givens_right = 0;
	Givens_left_tmp->deleteMat();	Givens_left_tmp = 0;
	D_tmp->deleteMat();				D_tmp = 0;
	return true;
}

/* Determination of a Matrix */
/* Return log-determination value of one square matrix */
/* Similar transformed  matrix is obtained via SVD for 
   matrix determination computation */

template<typename T>
float32 LogMatDet(Matrix<T> *A) {
	if (A==0) {
		cerr << "Det: The input Matrix should be initialized firstly" 
			 << endl;
		return -1;
	}
	Matrix<T> *U = 0; 
	Matrix<T> *V = 0;
	Matrix<T> *D = 0;
	svd(A, U, D, V);
	int16 d;
	d = (D->GetRowNum()>D->GetColNum()) ? D->GetColNum():D->GetRowNum();
	float32 tmp = 0.0;
	for (int i=0; i<d; i++) {
		tmp += log(D->data[i][i] + 1.0e-20);
	}
	U->deleteMat(); U = 0;
	D->deleteMat(); D = 0;
	V->deleteMat(); V = 0;
	return tmp;
}
}

#endif
