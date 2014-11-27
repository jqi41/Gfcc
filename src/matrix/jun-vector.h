#ifndef JUN_VECTOR_H_
#define JUN_VECTOR_H_	1

using namespace kaldi;

#include "../base/kaldi-common.h"

namespace jun {

template<typename T>
class Vector {
private:
	/* Size of the vector */
	int32 nSize;

	/* Resize the vector */
	void resize(int32 m);

public:
	/* Data of the vector */
	T *data;

	/* Classical creation */
	Vector() : data(0), nSize(0) {};

	/* Create a new vector */
	bool New(int32 n_dim);

	/* Compute the norm1 */
	T norm1(Vector<T> *weights = 0);

	/* Compute the norm2 */
	T norm2(Vector<T> * weights = 0);

	/* Compute the norm inf */
	T normInf();

	/* Delete the vector */
	void delVec();

	/* Output the element of vector */
	void print() const;

	/* Output the size of vector */
	int32 size() const;

	/* Output Data stream */
	T* OutDatStream();

	/* Zero the vector */
	void zero();

	virtual	~Vector();
};

template<typename T>
bool Vector<T>::New(int32 n_dim)
{
	nSize = n_dim;
	data = new T[n_dim];
	zero();
	return true;
}

template<typename T>
void Vector<T>::zero()
{
	if (data == NULL)
		return;
	else {
		for (int32 i=0; i<nSize; i++)
			data[i] = 0.0;
		return;
	}
}

template<typename T>
void Vector<T>::delVec()
{
	delete []data;
	data = 0;
	nSize = 0;
	return;
}

template<typename T>
void Vector<T>::resize(int32 m)
{
	T *_dat = new T[m];

	if (nSize < m)
		for (int32 i=0; i<nSize; i++)
			_dat[i] = data[i];
	else if (nSize > m)
		for (int32 i=0; i<m; i++)
			_dat[i] = data[i];

	delVec();
	data = new T[m];

	if (nSize < m) {
		for (int32 i=0; i<nSize; i++)
			data[i] = _dat[i];
		for (int32 i=nSize; i<m; i++)
			data[i] = 0.0;
	}
	else if (nSize > m) {
		for (int32 i=0; i<m; i++)
			data[i] = _dat[i];
	}
	nSize = m;
	delete _dat;
	return;
}

template<typename T>
inline T* Vector<T>::OutDatStream()
{
	return data;
}

template<typename T>
T Vector<T>::norm1(Vector<T> *weights)
{
	T sum = 0.0;
	T *_dat = data;
	if (weights) {
		T* ptr_w = weights->OutDatStream();
		for (int32 i=0; i<nSize; i++)
			sum += (*ptr_w++) * fabs(*_dat++);
	}
	else {
		for (int32 i=0; i<nSize; i++)
			sum += fabs(*_dat++);
	}
	delete []_dat;
	return sum;
}

template<typename T>
T Vector<T>::norm2(Vector<T> *weights)
{
	T sum = 0.0;
	T *_dat = data;
	if (weights) {
		T* ptr_w = weights->OutDatStream();
		for (int32 i=0; i<nSize; i++)
			sum += pow(((*ptr_w++) * fabs(*_dat++)), 2);
	}
	else {
		for (int32 i=0; i<nSize; i++)
			sum += pow(fabs(*_dat++), 2);
	}
	sum = sqrt(sum);

	return sum;
}

template<typename T>
T Vector<T>::normInf()
{
	T *_dat = data;
	T max_val = fabs(*_dat++);
	for (int32 i=1; i<nSize; i++)
	{
		T z = fabs(*_dat++);
		if (max_val < z)
			max_val = z;
	}
	return max_val;
}

template<typename T>
void Vector<T>::print() const
{
	for (int32 i=0; i<nSize; i++)
		std::cout << data[i] << " ";
	return;
}

template<typename T>
int32 Vector<T>::size() const
{
	return nSize;
}

template<typename T>
Vector<T>::~Vector() {
	delete []data;
	nSize = 0;
}
}
#endif
