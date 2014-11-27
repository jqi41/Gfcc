#ifndef JUN_FEATGFCC_H_
#define JUN_FEATGFCC_H_	1

#include "../base/kaldi-common.h"
#include "../util/parse-options.h"
#include "../feat/jun-featbase.h"
#include "../feat/jun-wave.h"
#include "../matrix/jun-matrix.h"
#include "../matrix/jun-vector.h"

using namespace kaldi;
using namespace std;

namespace jun {

#define MAXGLOSLEN 2048
#define DB	60.0
#define erb(x) ( 24.7 * (4.37e-3 * (x) + 1.0) )
#define VERY_SMALL_VALUE	1e-200
#define BW_CORRECTION	1.0190

struct Gfcc_Options {
public:
    BaseFloat sampPeriod;
    int32 StaticDimNum;
    int32 frmSize;
    int32 Increment;
    int32 Smp_Freq;
    int32 Delta_Order;
	int32 nChans;
	int32 Low_cf;
	int32 High_cf;
    bool  HasC0;
	bool  WavHead;
    string FeatNormKind;
	string FeatFn;
	string WavFn;

    explicit Gfcc_Options() : StaticDimNum(12), frmSize(400), 
		Increment(160), Smp_Freq(16000), Delta_Order(1), nChans(128),
		HasC0(false), FeatNormKind("CMS"), WavHead(true),
		Low_cf(80), High_cf(5000) {} 

    void Register(kaldi::ParseOptions *po) {
        po->Register("FrmSize", &frmSize, "Frame window length");
        po->Register("Increment", &Increment, "Increment target rate");
        po->Register("Smp_Freq", &Smp_Freq, "Sampling rate");
		po->Register("Low_Freq", &Low_cf, "Low center frequency");
		po->Register("High_Freq", &High_cf, "High center frequency");
        po->Register("Delta_Order", &Delta_Order, "Delta order");
        po->Register("HasC0", &HasC0, "Has C0 feature");
		po->Register("NChans", &nChans, "Number of channels");
        po->Register("StaticDimNum", &StaticDimNum, "Static dimension number");
        po->Register("FeatNormKind", &FeatNormKind, "Feature normalized kind");
		po->Register("AddWavHead", &WavHead, "Had wave header");
    }
};

struct Gammatone {
	BaseFloat *cf;		/* Center frequency for Gammatone filter bank */
	BaseFloat *bf;		/* Bandwidth of a given center frequency */
	BaseFloat cf_start;	/* Start center frequency of gammatone filter */
	BaseFloat cf_end;	/* End center frequency of gammatone filter */
	int fs;			/* Frequency sampling rate */
	explicit Gammatone() : cf_start(80), cf_end(5000), fs(16000), 
		cf(0), bf(0)
	{}
	void Clear() {
		delete []cf; cf = 0;
		delete []bf; bf = 0;
	}
};

struct Gfcc_feat {
	Vector<BaseFloat> static_feat;
	Vector<BaseFloat> dynamic_feat;
	Vector<BaseFloat> accelerate_feat;
	BaseFloat	static_C0;
	BaseFloat	dynamic_C0;
	BaseFloat	accelerate_C0;
	explicit Gfcc_feat() : static_C0(0), dynamic_C0(0),
			accelerate_C0(0) {}
};

class feat_gfcc : public feat_base {
private:
	Audio *audio_fmt;
	Gammatone *_gt;
	Gfcc_feat *gfcc_feat;
	Matrix<BaseFloat>	*out_gt;
	Matrix<BaseFloat>	*out_cg;
	Vector<BaseFloat>	*out_c0;
	HTKhdr	*htk_head;
	Gfcc_Options gfcc_cfg;
	int32 nFrames;
public:
	feat_gfcc() : audio_fmt(0), _gt(0), gfcc_feat(0), out_c0(0),
		out_gt(0), out_cg(0), htk_head(0), nFrames(0) {};
	virtual ~feat_gfcc() {}
	void GammatoneInit();
	bool GammatoneFilter(Audio *_in);
	void NormFeat();
	bool Gfcc_Gen(Gfcc_Options &gfcc_cfg);
	void CalDynamicFeat();
	void CalAccelerateFeat();
	bool FBank2C0();
	bool WriteFtrDat(const string FeatFn, const string WavFn);
	void ObsPrint(int nfrm, bool isOut) const;
	void Clear();
};
}
#endif