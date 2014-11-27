#include "../feat/jun-featGfcc.h"

namespace jun {

void feat_gfcc::GammatoneInit() {
	if (_gt==0) {
		_gt = new Gammatone;
		_gt->bf = new BaseFloat[gfcc_cfg.nChans];
		_gt->cf = new BaseFloat[gfcc_cfg.nChans];
	}
	_gt->fs = gfcc_cfg.Smp_Freq;
	BaseFloat erb_lo = Frq2Erb(_gt->cf_start);
	BaseFloat erb_hi = Frq2Erb(_gt->cf_end);
	BaseFloat erb_aver = (erb_hi - erb_lo) / (gfcc_cfg.nChans - 1);
	BaseFloat tmp;
	int  i=0;
	while (i<gfcc_cfg.nChans) {
		tmp = i * erb_aver + erb_lo;
		_gt->cf[i] = Erb2Frq(tmp);
		_gt->bf[i] = 1.019 * (24.7 + 0.108*_gt->cf[i]);
		i++;
	}
	return;
}

bool feat_gfcc::GammatoneFilter(Audio *in) {
	if (_gt==0)
		GammatoneInit();
	if (in==0) {
		cout << "Input data should be initialized" << endl;
		return false;
	}
	int32 nSamples = in->dataLen / sizeof(int16);
	if (out_gt==0) {
		out_gt = new Matrix<BaseFloat>;
		out_gt->CreateMat(gfcc_cfg.nChans, nSamples);
	}
	Matrix<BaseFloat> *senv = new Matrix<BaseFloat>;
	senv->CreateMat(gfcc_cfg.nChans, nSamples);

	double a, tpt, tptbw;
	double p0r, p1r, p2r, p3r, p4r;
	double p0i, p1i, p2i, p3i, p4i;
	double a1, a2, a3, a4, a5, u0r, u0i;
	double qcos, qsin, oldcs, coscf, sincf;
//	double oldphase = 0.0;
	BaseFloat floor_value = 1e-5;
	BaseFloat gain;
	int nSamples_padded = gfcc_cfg.Increment * nFrames;
	if (out_cg==0) {
		out_cg = new Matrix<BaseFloat>;
		out_cg->CreateMat(gfcc_cfg.nChans, nFrames);
	}
	tpt = (M_PI + M_PI) / _gt->fs;

	for (int n=0; n<gfcc_cfg.nChans; n++) {
		tptbw = tpt * erb(_gt->cf[n]) * BW_CORRECTION;
		a = exp(-tptbw);
		gain = tptbw * tptbw * tptbw * tptbw / 3;

		/* Variables initialization for each different channel */
		a1 = 4.0*a;	a2 = -6.0*a*a;	a3 = 4.0*a*a*a;	a4 = -a*a*a*a;	a5 = 4.0*a*a;
		p0r = 0.0; p1r = 0.0; p2r = 0.0; p3r = 0.0; p4r = 0.0;
		p0i = 0.0; p1i = 0.0; p2i = 0.0; p3i = 0.0; p4i = 0.0;
		coscf = cos(tpt * _gt->cf[n]);
		sincf = sin(tpt * _gt->cf[n]);
		qcos = 1;	qsin = 0;

		for(int t=0; t<nSamples; t++) {
			/* Filter part 1 & shift down to d.c. */
			p0r = qcos*in->data[t] + a1*p1r + a2*p2r + a3*p3r + a4*p4r;
			p0i = qsin*in->data[t] + a1*p1i + a2*p2i + a3*p3i + a4*p4i;

			/* Clip coefficients to stop them from becoming too close to zero */
			if (std::fabs(p0r) < VERY_SMALL_VALUE)
				p0r = 0.0;
			if (std::fabs(p0i) < VERY_SMALL_VALUE)
				p0i = 0.0;

			/* Filter part 2 */
			u0r = p0r + a1*p1r + a5*p2r;
			u0i = p0i + a1*p1i + a5*p2i;

			/* Update filter results */
			p4r = p3r; p3r = p2r; p2r = p1r; p1r = p0r;
			p4i = p3i; p3i = p2i; p2i = p1i; p1i = p0i;

			/*==========================================*/
			/* Basilar membrane response
			/* 1/ shift up in frequency first: (u0r+i*u0i) * exp(i*tpt*cf*t) = (u0r+i*u0i) * (qcos + i*(-qsin))
			/* 2/ take the real part only: bm = real(exp(j*wcf*kT).*u) * gain;
			/*==========================================*/
			out_gt->data[n][t] = ( u0r * qcos + u0i * qsin  ) * gain;
			senv->data[n][t] = sqrt( u0r * u0r + u0i * u0i ) * gain;

			oldcs = qcos;
			qcos = coscf * oldcs + sincf * qsin;
			qsin = coscf * qsin - sincf * oldcs;
		}
		for (int t=nSamples; t<nSamples_padded; t++) {
			p0r = a1*p1r + a2*p2r + a3*p3r + a4*p4r;
			p0i = a1*p1i + a2*p2i + a3*p3i + a4*p4i;
			u0r = p0r + a1*p1r + a5*p2r;
			u0i = p0i + a1*p1i + a5*p2i;
			p4r = p3r; p3r = p2r; p2r = p1r; p1r = p0r;
			p4i = p3i; p3i = p2i; p2i = p1i; p1i = p0i;

			senv->data[n][t] = sqrt( u0r * u0r + u0i * u0i ) * gain;
		}
		for (int m=0; m<nFrames; m++) {
			for (int k=0; k<gfcc_cfg.frmSize; k++)
				out_cg->data[n][m] += senv->data[n][k+m*gfcc_cfg.Increment];
			out_cg->data[n][m] /= gfcc_cfg.frmSize;
			out_cg->data[n][m] = pow((BaseFloat)out_cg->data[n][m], (BaseFloat)0.3);
			if (out_cg->data[n][m] < floor_value)
				out_cg->data[n][m] = floor_value;
			out_cg->data[n][m] = log(out_cg->data[n][m]);
		}
	}
	senv->deleteMat();
	delete senv;
	return true;
}

void feat_gfcc::CalDynamicFeat() {
	BaseFloat sum = 0.0;
	BaseFloat sum_E = 0.0;
	int16 i=1;
	int16 div = 0;
	int16 diff_back;
	int16 diff_front;
	BaseFloat back, head;
	BaseFloat back_E, head_E;
	int winDur = 2;
	for (int n=1; n<=nFrames; n++) {
		for (int m=0; m<gfcc_cfg.StaticDimNum; m++) {
			while (i<=winDur) {
				diff_back = n - i - 1;
				diff_front = n + i - 1;
				if (diff_back>0) {
					back = gfcc_feat[diff_back].static_feat.data[m];
					back_E = gfcc_feat[diff_back].static_C0;
				}
				else {
					back = gfcc_feat[0].static_feat.data[m];
					back_E = gfcc_feat[0].static_C0;
				}
				if (diff_front<nFrames) {
					head = gfcc_feat[diff_front].static_feat.data[m];
					head_E = gfcc_feat[diff_front].static_C0;
				}
				else {
					head = gfcc_feat[nFrames-1].static_feat.data[m];
					head_E = gfcc_feat[nFrames-1].static_C0;
				}
				sum += i*(head - back);
				sum_E += i*(head_E - back_E);
				div = div + i*i;
				i++;
			}
			sum /= 2*div;
			sum_E /= 2*div;
			gfcc_feat[n-1].dynamic_feat.data[m] = sum;
			gfcc_feat[n-1].dynamic_C0 = sum_E;
			sum = 0.0;	sum_E = 0.0;
			i = 1;
			div = 0;
		}
	}
	return;
}

void feat_gfcc::CalAccelerateFeat() {
	BaseFloat sum = 0.0;
	BaseFloat sum_E = 0.0;
	int i=1;
	int div = 0;
	int diff_back;
	int diff_front;
	BaseFloat back, head;
	BaseFloat back_E, head_E;
	int winDur = 2;
	for (int n=1; n<=nFrames; n++) {
		for (int m=0; m<gfcc_cfg.StaticDimNum; m++) {
			while (i<=winDur) {
				diff_back = n - i - 1;
				diff_front = n + i - 1;
				if (diff_back>0) {
					back = gfcc_feat[diff_back].dynamic_feat.data[m];
					back_E = gfcc_feat[diff_back].dynamic_C0;
				}
				else {
					back = gfcc_feat[0].dynamic_feat.data[m];
					back_E = gfcc_feat[0].dynamic_C0;
				}
				if (diff_front<nFrames) {
					head = gfcc_feat[diff_front].dynamic_feat.data[m];
					head_E = gfcc_feat[diff_front].dynamic_C0;
				}
				else {
					head = gfcc_feat[nFrames-1].dynamic_feat.data[m];
					head_E = gfcc_feat[nFrames-1].dynamic_C0;
				}
				sum += i*(head - back);
				sum_E += i*(head_E - back_E);
				div = div + i*i;
				i++;
			}
			sum /= 2*div;
			sum_E /= 2*div;
			gfcc_feat[n-1].accelerate_feat.data[m] = sum;
			gfcc_feat[n-1].accelerate_C0 = sum_E;
			sum = 0.0;	sum_E = 0.0;
			i = 1;
			div = 0;
		}
	}
	return;
}

void feat_gfcc::NormFeat()
{
	BaseFloat p=0.0;
	BaseFloat mean=0.0;
	BaseFloat var=0.0;
	if (!gfcc_cfg.FeatNormKind.compare("CMS")) {
		for(int i=0; i<gfcc_cfg.StaticDimNum; i++) {
			for(int k=0; k<nFrames; k++) {
				mean += gfcc_feat[k].static_feat.data[i];
			}
			mean /= nFrames;
			for(int k=0; k<nFrames; k++) {
				gfcc_feat[k].static_feat.data[i] -= mean;
			}
			mean = 0.0;
		}
	}
	if (!gfcc_cfg.FeatNormKind.compare("MVN")) {
		for(int i=0; i<gfcc_cfg.StaticDimNum; i++) {
			for(int k=0; k<nFrames; k++) {
				p += (gfcc_feat[k].static_feat.data[i])*(gfcc_feat[k].static_feat.data[i]);
				mean += gfcc_feat[k].static_feat.data[i];
			}
			p /= (nFrames-1);
			mean /= nFrames;
			var = sqrt(p-mean*mean);
			for(int k=0; k<nFrames; k++)
				gfcc_feat[k].static_feat.data[i] = (gfcc_feat[k].static_feat.data[i] - mean)/var;
			p = mean = var = 0.0;
		}
	}
	return;
}

bool feat_gfcc::WriteFtrDat(const std::string FeatFn, const std::string WavFn) {
	BaseFloat tmp;
	std::ofstream *out = new std::ofstream;
	out->open(FeatFn.c_str(), std::ios::trunc|std::ios::out|std::ios::binary);

	if (out->fail()) {
		std::cout << "The output stream could be not created !" << std::endl;
		return false;
	}

	if (htk_head==0)
		htk_head = new HTKhdr;
	out->seekp(0, std::ios::beg);
	int dim = gfcc_cfg.StaticDimNum;

	if (gfcc_cfg.HasC0)
		dim++;

	if (gfcc_cfg.Delta_Order==1)
		dim *= 2;

	else if (gfcc_cfg.Delta_Order==2)
		dim *= 3;

	htk_head->nSamples = nFrames;
	htk_head->sampPeriod = gfcc_cfg.Increment * (1e7 / gfcc_cfg.Smp_Freq);
	htk_head->parmName.append("USER");
	htk_head->sampSize = sizeof(BaseFloat) * dim;
	WriteHTKHeader(*out, htk_head, true);

	for (int i=0; i<nFrames; i++) {
		for (int k=0; k<gfcc_cfg.StaticDimNum; k++) {
			tmp = gfcc_feat[i].static_feat.data[k];
			Swap32((int32*)&tmp);
			out->write((char*)&tmp, sizeof(tmp));
		}
		if (gfcc_cfg.HasC0) {
			tmp = gfcc_feat[i].static_C0;
			Swap32((int32*)&tmp);
			out->write((char*)&tmp, sizeof(tmp));
		}
		if (gfcc_cfg.Delta_Order > 0) {
			for (int k=0; k<gfcc_cfg.StaticDimNum; k++) {
				tmp = gfcc_feat[i].dynamic_feat.data[k];
				Swap32((int32*)&tmp);
				out->write((char*)&tmp, sizeof(tmp));
			}
			if (gfcc_cfg.HasC0) {
				tmp = gfcc_feat[i].dynamic_C0;
				Swap32((int32*)&tmp);
				out->write((char*)&tmp, sizeof(tmp));
			}
			if (gfcc_cfg.Delta_Order==2) {
				for (int k=0; k<gfcc_cfg.StaticDimNum; k++) {
					tmp = gfcc_feat[i].accelerate_feat.data[k];
					Swap32((int32*)&tmp);
					out->write((char*)&tmp, sizeof(tmp));
				}
				if (gfcc_cfg.HasC0) {
					tmp = gfcc_feat[i].accelerate_C0;
					Swap32((int32*)&tmp);
					out->write((char*)&tmp, sizeof(tmp));
				}
			}
		}
	}
	std::cout << "Feature extraction: " << WavFn.c_str() << " --> " <<
		FeatFn.c_str() << std::endl;
	out->close();
	delete out;
	return true;
}

bool feat_gfcc::FBank2C0() {
	if (out_cg==0) {
		cerr << "Input cochleagram should be accessed at first." << endl;
		return false;
	}
	if (out_c0==0) {
		out_c0 = new Vector<BaseFloat>;
		out_c0->New(out_cg->GetColNum());
	}
	BaseFloat mfnorm = sqrt(2.0/(BaseFloat)gfcc_cfg.nChans);
	BaseFloat sum = 0.0;
	for (int i=0; i<out_cg->GetColNum(); i++) {
		for (int j=0; j<gfcc_cfg.nChans; j++) {
			sum += out_cg->data[j][i];
		}
		out_c0->data[i] = sum * mfnorm;
		sum = 0.0;
	}
	return true;
}

void feat_gfcc::ObsPrint(const int m, bool OutSmp) const {
	if (OutSmp) {
		std::cout << "Static features:" << std::endl;
		for (int i=0; i<gfcc_cfg.StaticDimNum; i++)
			std::cout << gfcc_feat[m].static_feat.data[i] << std::endl;
		if (gfcc_cfg.HasC0)
			std::cout << gfcc_feat[m].static_C0 << std::endl;
		if (gfcc_cfg.Delta_Order==2) {
			std::cout << "\nDynamic features:" << std::endl;
			for (int i=0; i<gfcc_cfg.StaticDimNum; i++)
				std::cout << gfcc_feat[m].dynamic_feat.data[i] << std::endl;
			if (gfcc_cfg.HasC0)
				std::cout << gfcc_feat[m].dynamic_C0 << std::endl;
			std::cout << "\nAccelerate features" << std::endl;
			for (int i=0; i<gfcc_cfg.StaticDimNum; i++)
				std::cout << gfcc_feat[m].accelerate_feat.data[i] << std::endl;
			if (gfcc_cfg.HasC0)
				std::cout << gfcc_feat[m].accelerate_C0 << std::endl;
		}
		else if (gfcc_cfg.Delta_Order==1) {
			std::cout << "\nDynamic features:" << std::endl;
			for (int i=0; i<gfcc_cfg.StaticDimNum; i++)
				std::cout << gfcc_feat[m].dynamic_feat.data[i] << std::endl;
			if (gfcc_cfg.HasC0)
				std::cout << gfcc_feat[m].dynamic_C0 << std::endl;
		}
	}
	return;
}

void feat_gfcc::Clear() {
	if (audio_fmt!=0) {
		delete []audio_fmt->data;
		audio_fmt = 0;
	}
	if (out_c0!=0) {
		out_c0->delVec();	out_c0 = 0;
	}
	if (_gt!=0) { _gt->Clear();	_gt = 0; }
	if (gfcc_feat!=0) {
		for (int n=0; n<nFrames; n++) {
			gfcc_feat[n].static_feat.delVec();
			if (gfcc_cfg.Delta_Order==1)
				gfcc_feat[n].dynamic_feat.delVec();
			else if (gfcc_cfg.Delta_Order==2) {
				gfcc_feat[n].dynamic_feat.delVec();
				gfcc_feat[n].accelerate_feat.delVec();
			}
		}
		delete []gfcc_feat;
		gfcc_feat = 0;
	}
	if (out_gt!=0) { out_gt->deleteMat();	out_gt=0; }
	if (out_cg!=0) { out_cg->deleteMat();	out_cg=0; }
	if (htk_head!=0) { delete htk_head; htk_head=0; }
	gfcc_cfg.FeatNormKind.clear();
	gfcc_cfg.FeatFn.clear();
	gfcc_cfg.WavFn.clear();
}

bool feat_gfcc::Gfcc_Gen(Gfcc_Options &_gfcc_cfg) {
	gfcc_cfg = _gfcc_cfg;
	Wave src(gfcc_cfg.WavFn.c_str());
	if (audio_fmt==0)
		audio_fmt = new Audio;
	src.ReadWavData(*audio_fmt, gfcc_cfg.WavHead);
	nFrames = static_cast<int32>(((audio_fmt->dataLen/sizeof(short)) \
		- gfcc_cfg.frmSize) / gfcc_cfg.Increment) + 1;
	if (gfcc_feat==0)
		gfcc_feat = new Gfcc_feat[nFrames];
	GammatoneFilter(audio_fmt);
	FBank2C0();

	/* DCT Computation */
	BaseFloat tmp;
	BaseFloat fCoef = sqrt(2.0 / gfcc_cfg.nChans);
	for (int n=0; n<nFrames; n++) {
		gfcc_feat[n].static_feat.New(gfcc_cfg.StaticDimNum);
		for (int m=1; m<=gfcc_cfg.StaticDimNum; m++) {
			tmp = 0.0;
			for (int k=1; k<=gfcc_cfg.nChans; k++)
				tmp = tmp + cos(m*M_PI*(2*k-1)/(2*gfcc_cfg.nChans)) * out_cg->data[k-1][n];
			gfcc_feat[n].static_feat.data[m-1] = tmp;
			gfcc_feat[n].static_feat.data[m-1] *= fCoef;
		}
	}
	if (gfcc_cfg.HasC0) {
		for (int i=0; i<nFrames; i++) {
			gfcc_feat[i].static_C0 = out_c0->data[i];
		}
	}

	/* Compute dynamic and accelerate features */
	if (gfcc_cfg.Delta_Order==1) {			/* for the case of MFCC_D */
		for (int i=0; i<nFrames; i++) {
			gfcc_feat[i].dynamic_feat.New(gfcc_cfg.StaticDimNum);
		}
		CalDynamicFeat();
	}
	else if(gfcc_cfg.Delta_Order==2) {		/* for the case of MFCC_D_A */
		for (int i=0; i<nFrames; i++) {
			gfcc_feat[i].dynamic_feat.New(gfcc_cfg.StaticDimNum);
			gfcc_feat[i].accelerate_feat.New(gfcc_cfg.StaticDimNum);
		}
		CalDynamicFeat();
		CalAccelerateFeat();
	}

	/* Whether do feature normalization */
	if ( gfcc_cfg.FeatNormKind.compare("CMS") ||
		 gfcc_cfg.FeatNormKind.compare("MVN") )
		NormFeat();

	WriteFtrDat(gfcc_cfg.FeatFn.c_str(), gfcc_cfg.WavFn.c_str());
	ObsPrint(0, false);
	Clear();
	return true;
}
}
