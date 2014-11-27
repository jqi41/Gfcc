#ifndef _JUN_STATS_INF_H
#define _JUN_STATS_INF_H

#include "../base/kaldi-common.h"
#include "../matrix/jun-math.h"

namespace jun {
struct Mean_set {
	Vector<float32> *mean_val;
	Vector<float32> *sum_val;
	int32 count;
	explicit Mean_set() : count(0), mean_val(0), sum_val(0) {}
	void Init(int32 _count) {
		if (mean_val==0)
			mean_val = new Vector<float32>;
		if (sum_val==0)
			sum_val = new Vector<float32>;
		count = _count;
		mean_val->New(count);
		sum_val->New(count);
	}
	void Delete() {	
		if (mean_val!=0) 
			mean_val->delVec();
		if (sum_val!=0) 
			sum_val->delVec();
		count = 0;
	}
};

struct Cov_set {
	Matrix<float32> *cov_val;
	Matrix<float32> *sum_val;
	int32 count;
	explicit Cov_set() : count(0), cov_val(0), sum_val(0) {}
	void Init(int32 _count) {
		count = _count;
		if (cov_val==0)
			cov_val = new Matrix<float32>;
		if (sum_val==0) 
			sum_val = new Matrix<float32>;
		cov_val->CreateMat(count, count);
		sum_val->CreateMat(count, count);
	}
	void Delete() {
		if (cov_val!=0) 
			cov_val->deleteMat();
		if (sum_val!=0)
			sum_val->deleteMat();
		count = 0;
	}
};

void Mean_val(Matrix<float32> *_mat, Mean_set *_mean) {
	float32 tmp = 0.0;
	int32 count = _mat->GetColNum();
	if (_mean==0) {
		_mean = new Mean_set;
		_mean->Init(_mat->GetRowNum());
	}
	for (int32 i=0; i<_mat->GetRowNum(); i++) {
		for (int32 j=0; j<_mat->GetColNum(); j++) {
			tmp += _mat->data[i][j];
		}
		_mean->sum_val->data[i] = tmp;
		_mean->mean_val->data[i] = _mean->sum_val->data[i] / count;
		tmp = 0.0;
	}
	return;
}

void Cov_vval(Matrix<float32> *_mat, Mean_set *_mean, Cov_set *_cov) {
	Vector<float32> tmp_b;
	Vector<float32> tmp_m;
	float32 res;
	if (_cov==0) {
		_cov = new Cov_set;
		_cov->Init(_mat->GetRowNum());
	}
	tmp_b.New(_mat->GetRowNum());
	tmp_m.New(_mat->GetRowNum());
	for (int32 n=0; n<_mat->GetRowNum(); n++) {
		for (int32 m=0; m<_mat->GetColNum(); m++) {
			tmp_b.data[m] = _mat->data[n][m] - _mean->mean_val->data[n];
		}
		for (int32 m=n; m<_mat->GetRowNum(); m++) {
			for (int32 j=0; j<_mat->GetColNum(); j++) {
				tmp_m.data[j] = _mat->data[m][j] - _mean->mean_val->data[m];
			}
			mxInProd<float32>(&tmp_b, &tmp_m, res);
			_cov->sum_val->data[n][m] = res;
			_cov->sum_val->data[m][n] = _cov->sum_val->data[n][m];
			_cov->cov_val->data[n][m] = res / (_mat->GetColNum()-1);
			_cov->cov_val->data[m][n] = _cov->cov_val->data[n][m];
			tmp_m.zero();
		}
		tmp_b.zero();
	}
	return;
}
}

#endif