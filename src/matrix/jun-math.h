#ifndef JUN_MATH_H_
#define JUN_MATH_H_ 1

#include "../base/kaldi-common.h"
#include "../matrix/jun-matrix.h"
#include "../matrix/jun-vector.h"

using namespace std;
using namespace kaldi;

namespace jun {

inline void Swap16(int16 *Short) {
	*Short  = ((*Short & 0x00ff) << 8) | ((*Short & 0xff00) >> 8);
}

inline void Swap32(int32 *Long) {
	*Long = ((*Long&0x000000ffL)<<24 )| ((*Long&0x0000ff00L)<<8 )| \
			((*Long&0x00ff0000L)>>8 ) | ((*Long&0xff000000L)>>24 );
}

inline int32 minab(int a, int b) {
	if (a > b)
		return b;
	else
		return a;
}

inline int32 sgn(int32 x) {
	if (x >= 0)
		return 1;
	else
		return -1;
}

/* Copy a matrix from another matrix */
template<typename T>
bool CopyMatrix(Matrix<T> *&mat, Matrix<T> *mat_tmp) {
	if (mat_tmp==0) {
		cout << "Temporary matrix should be initialzied." << endl;
		return false;
	}
	int32 col = mat_tmp->GetColNum();
	int32 row = mat_tmp->GetRowNum();

	if (mat==0) {
		mat = new Matrix<T>;
		mat->CreateMat(row, col);
	}
	for (int n=0; n<row; n++) {
		for (int m=0; m<col; m++) {
			mat->data[n][m] = mat_tmp->data[n][m];
		}
	}
	return true;
}

/* Copy a vector from a specified column of the Matrix */
template<typename T>
bool CopyVector(Matrix<T> *_mat, int col, Vector<T> *&_vec) {
	if (_mat==0) {
		cout << "The matrix should be initialized at first" << endl;
		return false;
	}
	if (_vec==0)
		_vec->New(_mat->GetRowNum());
	int size = _mat->GetRowNum();
	for (int i=0; i<size; i++) {
		_vec->data[i] = _mat->data[i][col];
	}
	return true;
}

/* Copy a row vector from a specified row of the Matrix */
template<typename T>
bool CopyRow(Matrix<T> *_mat, int row, Vector<T> *&_row) {
	if (_mat==0) {
		cout << "The matrix should be initialized at first" << endl;
		return false;
	}
	if (_row==0)
		_row->New(_mat->GetColNum());
	for (int i=0; i<_mat->GetColNum(); i++) {
		_row->data[i] = _mat->data[row][i];
	}
	return true;
}

/* Matrix addition */
template<typename T>
bool mxMatAddMat(Matrix<T> *mat1, Matrix<T> *mat2, Matrix<T> *&out) {
	if (mat1->GetColNum() != mat2->GetColNum() &&
		mat2->GetRowNum() != mat2->GetRowNum()) {
		cout << "The column and row number should be the same." << endl;
		return false;
	}
	if (out==0) {
		out = new Matrix<T>;
		out->CreateMat(mat1->GetRowNum(), mat1->GetColNum());
	}
	for (int n=0; n<mat1->GetRowNum(); n++) {
		for (int m=0; m<mat1->GetColNum(); m++) {
			out->data[n][m] = mat1->data[n][m] + mat2->data[n][m];
		}
	}
	return true;
}

/* Matrix subtraction */
template<typename T>
bool mxMatSubMat(Matrix<T>* mat1, Matrix<T>* mat2, Matrix<T>* &out) {
	if (mat1->GetColNum() != mat2->GetColNum() &&
		mat2->GetRowNum() != mat2->GetRowNum()) {
		cout << "The column and row number should be the same." << endl;
		return false;
	}
	if (out==0) {
		out = new Matrix<T>;
	}
	out->CreateMat(mat1->GetRowNum(), mat1->GetColNum());

	for (int n=0; n<mat1->GetRowNum(); n++) {
		for (int m=0; m<mat1->GetColNum(); m++) {
			out->data[n][m] = mat1->data[n][m] - mat2->data[n][m];
		}
	}
	return true;
}

/* Matrix multiplication */
template<typename T>
bool mxMatMulMat(Matrix<T> *mat1, Matrix<T> *mat2, Matrix<T> *&out) {
	if (mat1->GetColNum() != mat2->GetRowNum()) {
		cout << "The column and row number should be the same." << endl;
		return false;
	}
	if (out==0) {
		out = new Matrix<T>;
	}
	out->CreateMat(mat1->GetRowNum(), mat2->GetColNum());
	float32 tmp = 0.0;
	for (int n=0; n<mat1->GetRowNum(); n++) {
		for (int l=0; l<mat2->GetColNum(); l++) {
			for (int m=0; m<mat1->GetColNum(); m++) {
				tmp += mat1->data[n][m] * mat2->data[m][l];
			}
			out->data[n][l] = tmp;
			tmp = 0.0;
		}
	}
	return true;
}

/* Matrix transposition */
template<typename T>
void mxTrMat(Matrix<T> *in, Matrix<T> *&out) {
	if (out==0) {
		out = new Matrix<T>;
		out->CreateMat(in->GetColNum(), in->GetRowNum());
	}
	if (in != out) {
		for (int i=0; i<in->GetRowNum(); i++) {
			for (int j=0; j<in->GetColNum(); j++) {
				out->data[j][i] = in->data[i][j];
			}
		}
	}
	else {
		for (int i=1; i<in->GetRowNum(); i++) {
			for (int j=0; j<i; j++) {
				float32 tmp = in->data[i][j];
				in->data[i][j] = in->data[j][i];
				in->data[j][i] = tmp;
			}
		}
	}
	return;
}

/* Multiplication of a matrix and a vector */
template<typename T>
bool mxMatMulVec(Matrix<T> *mat, Vector<T> *vec, Vector<T> *&out) {
	if (mat->GetColNum() != vec->size()) {
		cout << "The column number should be the same as the row" << endl;
		return false;
	}
	float32 tmp = 0.0;
	for (int m=0; m<mat->GetRowNum(); m++) {
		for (int n=0; n<vec->size(); n++) {
			tmp += mat->data[m][n] * vec->data[n];
		}
		out->data[m] = tmp;
		tmp = 0.0;
	}
	return true;
}

/* Multiplication of a vector and a matrix */
template<typename T>
bool mxVecMulMat(Vector<T> *vec, Matrix<T> *mat, Vector<T> *&out) {
	if (vec->size() != mat->GetRowNum()) {
		cout << "The vector size should be the same as the row number of a matrix.";
		return false;
	}
	if (out==0) {
		out = new Vector<T>;
		out->New(vec->size());
	}
	float32 tmp = 0.0;
	for (int m=0; m<mat->GetColNum(); m++) {
		for (int n=0; n<vec->size(); n++) {
			tmp += vec->data[n] * mat->data[n][m];
		}
		out->data[m] = tmp;
		tmp = 0.0;
	}
	return true;
}

/* return a value of the quadratic equation d'* Cov * t */
template<typename T>
bool mxQuadValue(Vector<T> *d, Matrix<T> *Cov, Vector<T> *t, float32 &res) {
	Vector<T> *dd = 0;
	mxVecMulMat(d, Cov, dd);
	mxInProd(dd, t, res);
	return true;
}

/* Multiplicaiton of transposed matrix1 and matrix2 */
template<typename T>
bool mxTrMatMulMat(Matrix<T> *mat1, Matrix<T> *mat2, Matrix<T> *&out) {
	if (mat1->GetRowNum() != mat2->GetRowNum()) {
		cout << "The row of transposed matrix 1 should be the same as the row of matrix 2."
			<< endl;
		return false;
	}
	if (out==0)
		out = new Matrix<T>;
	out->CreateMat(mat1->GetColNum(), mat2->GetColNum());

	float32 tmp = 0.0;
	for (int n=0; n<mat1->GetColNum(); n++) {
		for (int l=0; l<mat2->GetColNum(); l++) {
			for (int m=0; m<mat1->GetRowNum(); m++) {
				tmp += mat1->data[n][m] * mat2->data[m][l];
			}
			out->data[n][l] = tmp;
			tmp = 0.0;
		}
	}
	return true;
}

/* Multiplication of Matrix1 and transposed matrix2 */
template<typename T>
bool mxMatMulTrMat(Matrix<T> *mat1, Matrix<T> *mat2, Matrix<T> *&out) {
	if (mat1->GetColNum() != mat2->GetColNum()) {
		cout << "The columns of mat1 should be the same of transposed mat2" << endl;
		return false;
	}
	if (out==0)
		out = new Matrix<T>;
	out->CreateMat(mat1->GetRowNum(), mat2->GetRowNum());

	float32 tmp = 0.0;
	for (int n=0; n<mat1->GetRowNum(); n++) {
		for (int l=0; l<mat2->GetRowNum(); l++) {
			for (int m=0; m<mat1->GetColNum(); m++) {
				tmp += mat1->data[n][m] * mat2->data[m][l];
			}
			out->data[n][l] = tmp;
			tmp = 0.0;
		}
	}
	return true;
}

/* Inner product of two vectors */
template<typename T>
bool mxInProd(Vector<T> *in1, Vector<T> *in2, T &out) {
	if (in1->size() != in2->size()) {
		cout << "The size of two vectors should be the same." << endl;
		return false;
	}
	float32 tmp = 0.0;
	for (int n=0; n<in1->size(); n++) {
		tmp += in1->data[n] * in2->data[n];
	}
	out = tmp;
	return true;
}

/* Outer product of two vectors */
template<typename T>
void mxOutProd(Vector<T> *in1, Vector<T> *in2, Matrix<T> *&mat) {
	if (mat==0)
		mat = new Matrix<T>;
	mat->CreateMat(in1->size(), in2->size());
	for (int n=0; n<in1->size(); n++) {
		for (int m=0; m<in2->size(); m++) {
			mat->data[n][m] = in1->data[n] * in2->data[m];
		}
	}
	return;
}

/* Vector addition */
template<typename T>
bool VectorAdd(Vector<T> *in1, Vector<T> *in2, Vector<T> *&out) {
	if (in1->size() != in2->size()) {
		cout << "Two vectors should be the same." << endl;
		return false;
	}
	if (out==0)
		out = new Vector<T>;
	out->New(in1->size());
	for (int i=0; i<in1->size(); i++)
		out->data[i] = in1->data[i] + in2->data[i];

	return true;
}

/* Vector subtraction */
template<typename T>
bool VectorSub(Vector<T> *in1, Vector<T> *in2, Vector<T> *&out) {
	if (in1->size() != in2->size()) {
		cout << "Two vectors should be the same." << endl;
		return false;
	}
	if (out==0) {
		out = new Vector<T>;
		out->New(in1->size());
	}
	for (int i=0; i<in1->size(); i++)
		out->data[i] = in1->data[i] - in2->data[i];

	return true;
}

/* TranVec2Mat */
template<typename T>
void TranVec2Mat(Vector<T> *vec, Matrix<T> *&mat) {
	if (mat==0)
		mat = new Matrix<T>;
	mat->CreateMat(vec->size(), vec->size());
	for (int n=0; n<vec->size(); n++) {
		mat->data[n][n] = vec->data;
	}
	return;
}

/* Multiplication of a constant scalar and matrix */
template<typename T>
bool ConstMulMat(T _const, Matrix<T> *&mat) {
	if (mat==0) {
		cout << "ConstMulMat: The matrix should be initialized at first"
			<< endl;
		return false;
	}
	for (int m=0; m<mat->GetRowNum(); m++) {
		for (int n=0; n<mat->GetColNum(); n++) {
			mat->data[m][n] *= _const;
		}
	}
	return true;
}

/* Matrix iterative multiplication */
template<typename T>
void mxLItrMatMul(Matrix<T> *&mat, Matrix<T> *&H_mat) {
	Matrix<T> *mat_tmp = 0;
	mxMatMulMat<T>(H_mat, mat, mat_tmp);
	mat->deleteMat();
	mat->CreateMat(mat_tmp->GetRowNum(), mat_tmp->GetColNum());
	for (int n=0; n<mat_tmp->GetRowNum(); n++)
		for (int m=0; m<mat_tmp->GetColNum(); m++)
			mat->data[n][m] = mat_tmp->data[n][m];
	mat_tmp->deleteMat();
	delete mat_tmp; mat_tmp = 0;
}

/* Matrix iterative multiplication */
template<typename T>
void mxRItrMatMul(Matrix<T> *&mat, Matrix<T> *&H_mat) {
	Matrix<T> *mat_tmp = 0;
	mxMatMulMat<T>(mat, H_mat, mat_tmp);
	mat->deleteMat();
	mat->CreateMat(mat_tmp->GetRowNum(), mat_tmp->GetColNum());
	for (int n=0; n<mat_tmp->GetRowNum(); n++)
		for (int m=0; m<mat_tmp->GetColNum(); m++)
			mat->data[n][m] = mat_tmp->data[n][m];
	mat_tmp->deleteMat();
	delete mat_tmp; mat_tmp = 0;
}

template<typename T>
bool CopySubMat(Matrix<T> *&mat, Matrix<T> *mat_tmp, bool eye=true) {
	if (mat_tmp==0) {
		cout << "The copy matrix should be initialized at first." << endl;
		return false;
	}
	if (mat==0) {
		mat = new Matrix<T>;
		mat->CreateMat(mat_tmp->GetRowNum(), mat_tmp->GetColNum());
	}
	if (mat->GetRowNum() < mat_tmp->GetRowNum() && mat->GetColNum() < mat_tmp->GetColNum()) {
		for (int i=0; i<mat->GetRowNum(); i++)
			for (int j=0; j<mat->GetColNum(); j++)
				mat->data[i][j] = mat_tmp->data[i][j];
	}
	else if (mat->GetRowNum() >= mat_tmp->GetRowNum() &&
		mat->GetColNum() >= mat_tmp->GetRowNum()) {
			for (int i=0; i<mat_tmp->GetRowNum(); i++)
				for (int j=0; j<mat_tmp->GetColNum(); j++)
					mat->data[i][j] = mat_tmp->data[i][j];
			float32 sub_row = mat->GetRowNum() - mat_tmp->GetRowNum();
			float32 sub_col = mat->GetColNum() - mat_tmp->GetColNum();
			if (sub_row==sub_col && eye) {
				for (int n=mat_tmp->GetColNum(); n<mat->GetColNum(); n++)
					mat->data[n][n] = 1.0;
			}
	}
	return true;
}

template<typename T>
void MeanVec(Matrix<T> *_mat, Vector<T> *&_mean) {
	if (_mat==0) {
		cerr << "The input matrix should be initialized firstly."
			<< endl;
		return;
	}
	float32 tmp = 0.0;
	if (_mean==0) {
		_mean = new Vector<float32>;
		_mean->New(_mat->GetRowNum());
	}
	float32 sum_tmp = 0.0;
	for (int32 i=0; i<_mat->GetRowNum(); i++) {
		for (int32 j=0; j<_mat->GetColNum(); j++)
			tmp = tmp + _mat->data[i][j];
		_mean->data[i] = tmp / _mat->GetColNum();
		tmp = 0.0;
	}
	return;
}

}

#endif
