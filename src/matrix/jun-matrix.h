#ifndef JUN_MATRIX_H_
#define JUN_MATRIX_H_	1

#include "../base/kaldi-common.h"

using namespace kaldi;
/*                                                                             */
/* ------------------------------ Class Package ------------------------------ */
/*																			   */
											                   
namespace jun {
template <typename T>
class Matrix {
	unsigned int row, col;					// row and column numbers
	void resize(int32 _row, int32 _col);	// Resize the matrix
public:
	/* Data stored in Matrix */
	T **data;

	/* Default creator */
	Matrix() : data(0), row(0), col(0) {}

	/* Create a matrix */
	Matrix(int32 _row, int32 _col);

	/* Create a matrix from data stream with one dimension */
	Matrix(T *_dat, int32 _row, int32 _col);

	/* Create a matrix from data stream with two dimensions */
	Matrix(T **_dat, int32 _row, int32 _col);

	virtual ~Matrix() {};

	/* Create Matrix from a calling function */
	void CreateMat(int32 _row, int32 _col);

	/* Create eye's matrix */
	void CreateEyeMat(int32 _dim);

	// Return data stream of row n
	T *GetRowStream(int32 n);

	// Return data stream of column m
	T *GetColStream(int32 m);

	/* Copy the matrix #mat# */
	void copy(Matrix& _mat);

	/* Zero the matrix */
	void zero();

	/* Get the row number of matrix */
	inline int GetRowNum();

	/* Get the column number of matrix */
	inline int GetColNum();

	/* Delete the existing Matrix */
	void deleteMat();

	/* Output Matrix elements */
	void print() const;

	/* Obtain a identity  */
	void InitIdentity();

};

/*                                                                           */
/* -------------------- Implementation of Class Matrix --------------------- */
/*															                 */

template<typename T>
Matrix<T>::Matrix(int32 _row, int32 _col) : row(_row), col(_col)
{
	data = new T*[row];
	for (int32 i=0; i<row; i++)
		data[i] = new T[col];
}

template<typename T>
Matrix<T>::Matrix(T *_dat, int32 _row, int32 _col) :
		row(_row), col(_col)
{
	int32 m = 0;
	data = new T*[row];
	for (int32 i=0; i<row; i++)
		data[i] = new T[col];

	for(int32 i=0; i<row; i++)
		for(int32 j=0; j<col; j++)
			data[i][j] = _dat[m++];
}

template<typename T>
Matrix<T>::Matrix(T **_dat, int32 _row, int32 _col) :
		row(_row), col(_col)
{
	data = new T*[row];
	for (int32 i=0; i<row; i++)
		data[i] = new T[col];

	for(int32 i=0; i<row; i++)
		for(int32 j=0; j<col; j++)
			data[i][j] = _dat[i][j];
}

template<typename T>
void Matrix<T>::zero()
{
	for(size_t i=0; i<row; i++)
		for(size_t j=0; j<col; j++)
			data[i][j] = 0.0;
	return;
}

template<typename T>
inline int Matrix<T>::GetRowNum()
{
	return row;
}

template<typename T>
inline int Matrix<T>::GetColNum()
{
	return col;
}

template<typename T>
T *Matrix<T>::GetRowStream(int32 n)
{
	T *dataStream = new T[col];
	for (int32 j=0; j<col; j++)
		dataStream[j] = data[n][j];
	return dataStream;
}

template<typename T>
T *Matrix<T>::GetColStream(int32 m)
{
	T *dataStream = new T[row];
	for (int32 j=0; j<row; j++)
		dataStream[j] = data[j][m];
	return dataStream;
}

template<typename T>
void Matrix<T>::copy(Matrix& _mat)
{
	zero();

	if ((row==_mat.GetRowNum())&&(col==_mat.GetColNum()))
	{
		for (int32 i=0; i<row; i++) {
			data[i] = _mat.GetRowStream(i);
		}
	}
	else {
		resize(_mat.GetRowNum(), _mat.GetColNum());
		for (int32 i=0; i<row; i++) {
			data[i] = _mat.GetRowStream(i);
		}
	}

	return;
}

template<typename T>
void Matrix<T>::print() const
{
	using std::cout;
	using std::endl;
	for(int32 i=0; i<row; i++) {
		for(size_t j=0; j<col; j++)
			cout << data[i][j] << " ";
		cout << endl;
	}
}

template<typename T>
void Matrix<T>::deleteMat()
{
	for(size_t i=0; i<row; i++)
		delete [] data[i];
	delete [] data;
	data = 0;
}

template<typename T>
void Matrix<T>::resize(int32 _row, int32 _col)
{
	deleteMat();
	row = _row;
	col = _col;
	data = new T*[row];
	for (int32 i=0; i<row; i++) {
		data[i] = new T[col];
	}
}

template<typename T>
void Matrix<T>::CreateMat(int32 _row, int32 _col)
{
	data = new T*[_row];
	for (int32 i=0; i<_row; i++)
		data[i] = new T[_col];
	row = _row;
	col = _col;
	zero();
	return;
}

template<typename T>
void Matrix<T>::InitIdentity() {
	int size = minab(GetRowNum(), GetColNum());
	for (int i=0; i<size; i++) {
		data[i][i] = 1.0;
	}
	return;	
}

template<typename T>
void Matrix<T>::CreateEyeMat(int32 _dim) {
	data = new T*[_dim];
	for (int32 i=0; i<_dim; i++)
		data[i] = new T[_dim];
	row = _dim;
	col = _dim;
	zero();
	for (int i=0; i<_dim; i++) 
		data[i][i] = 1.0;
	return;
}

}

#endif
