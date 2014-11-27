#ifndef JUN_FEATBASE_H_
#define JUN_FEATBASE_H_	1

#include "../base/kaldi-common.h"
#include "../base/jun-label.h"
#include "../matrix/jun-math.h"

using namespace kaldi;

namespace jun {

const int MAXSTRLEN = 1024;

struct FrmDat {
	BaseFloat	*_real;
	BaseFloat	*_img;
	BaseFloat	*_power;
    int32	frmSize;
	explicit FrmDat() : frmSize(0), _real(0), _img(0),
		_power(0) {}
	void Clear() {
		if (_real==0)	{ delete []_real; _real=0;}
		if (_img==0)	{ delete []_img; _img=0;}
		if (_power==0)	{ delete []_power; _power=0;}
		frmSize = 0;
	}
};

class feat_base {
protected:
    void GenHamWindow(BaseFloat *ham, int32 winSize);
	int16 GetParmKind(const std::string &szParmKind);
public:
	bool WriteHTKHeader(std::ofstream &out, HTKhdr *header, bool Swap);
    void rfft(FrmDat *_signal);
    BaseFloat Frq2Mel(const BaseFloat freq);
    BaseFloat Mel2Frq(const BaseFloat mel);
	BaseFloat Frq2Erb(const BaseFloat freq);
	BaseFloat Erb2Frq(const BaseFloat erb);
};

}

#endif
