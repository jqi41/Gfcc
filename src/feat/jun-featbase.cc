#include "../feat/jun-featbase.h"

using namespace kaldi;
using namespace std;

namespace jun {

void feat_base::GenHamWindow(BaseFloat *ham, int32 hamWinSize) {
    if (ham==0)
        ham = new BaseFloat[hamWinSize];

    BaseFloat a = 2 * M_PI / (hamWinSize-1);
    for (int i=0; i<hamWinSize; i++)
        ham[i] = 0.54 - 0.46 * cos(a*(i));
}

int16 feat_base::GetParmKind(const std::string &szParmKind)
{
	int16 sKind = 0;

	/* Base Type */
	if (NULL != strstr(szParmKind.c_str(), "WAVEFORM"))				/* sampled waveform */
		sKind |= 0x0;
	else if (NULL != strstr(szParmKind.c_str(), "LPC"))				/* linear prediction filter coefficients */
		sKind |= 0x1;
	else if (NULL != strstr(szParmKind.c_str(), "LPREFC"))			/* linear prediction reflection coefficients */
		sKind |= 0x2;
	else if (NULL != strstr(szParmKind.c_str(), "LPCEPSTRA"))		/* LPC cepstral coefficients */
		sKind |= 0x3;
	else if (NULL != strstr(szParmKind.c_str(), "LPDELCEP"))		/* LPC cepstra plus delta coefficients */
		sKind |= 0x4;
	else if (NULL != strstr(szParmKind.c_str(), "IREFC"))			/* LPC reflection coef in 16 bit integer format */
		sKind |= 0x5;
	else if (NULL != strstr(szParmKind.c_str(), "MFCC"))			/* mel-frequency cepstral coefficients */
		sKind |= 0x6;
	else if (NULL != strstr(szParmKind.c_str(), "FBANK"))			/* log mel-filter bank channel outputs */
		sKind |= 0x7;
	else if (NULL != strstr(szParmKind.c_str(), "MELSPEC"))			/* linear mel-filter bank channel outputs */
		sKind |= 0x8;
	else if (NULL != strstr(szParmKind.c_str(), "USER"))			/* user defined sample kind */
		sKind |= 0x9;
	else if (NULL != strstr(szParmKind.c_str(), "DISCRETE"))		/* vector quantised data */
		sKind |= 0xA;
	else if (NULL != strstr(szParmKind.c_str(), "PLP"))				/* PLP cepstral coefficients */
		sKind |= 0xB;
	else
		cout << "Unknown parameter kind" << endl;

	/* Extension */
	if (NULL != strstr(szParmKind.c_str(), "_E"))		/* has energy */
		sKind |= 000100;
	if (NULL != strstr(szParmKind.c_str(), "_N"))		/* absolute energy suppressed */
		sKind |= 000200;
	if (NULL != strstr(szParmKind.c_str(), "_D"))		/* has delta coefficients */
		sKind |= 000400;
	if (NULL != strstr(szParmKind.c_str(), "_A"))		/* has acceleration coefficients */
		sKind |= 001000;
	if (NULL != strstr(szParmKind.c_str(), "_C"))		/* is compressed */
		sKind |= 002000;
	if (NULL != strstr(szParmKind.c_str(), "_Z"))		/* has zero mean static coef. */
		sKind |= 004000;
	if (NULL != strstr(szParmKind.c_str(), "_K"))		/* has CRC chechsum */
		sKind |= 010000;
	if (NULL != strstr(szParmKind.c_str(), "_0"))		/* has 0'th cepstral coef. */
		sKind |= 020000;

	return sKind;
}

bool feat_base::WriteHTKHeader(std::ofstream &out, HTKhdr *header, bool Swap) {
	if (!out.is_open()) {
		std::cout << "No given pointer !" << std::endl;
		return false;
	}
	int16 sampKind = GetParmKind(header->parmName);
	int16 sampSize = header->sampSize;
	int32 nSamples = header->nSamples;
	int32 sampPeriod = header->sampPeriod;
	if ( Swap )
	{
		Swap32(&nSamples);
		Swap32(&sampPeriod);
		Swap16(&sampSize);
		Swap16(&sampKind);
	}
	out.write((char*)&nSamples, sizeof(int32));
	out.write((char*)&sampPeriod, sizeof(int32));
	out.write((char*)&sampSize, sizeof(int16));
	out.write((char*)&sampKind, sizeof(int16));
	out.seekp(0, std::ios::end);
	return true;
}

BaseFloat feat_base::Mel2Frq(const BaseFloat _mel) {
    return 700.0 * (exp(_mel/1127.0) - 1.0);
}

BaseFloat feat_base::Frq2Mel(const BaseFloat _freq) {
    return 1127.0 * log(1.0+_freq/700.0);
}

BaseFloat feat_base::Erb2Frq(const BaseFloat _erb) {
	return (exp(_erb/(21.4*log10(2.71828)))-1.0) * 1e3 / 4.37;
}

BaseFloat feat_base::Frq2Erb(const BaseFloat _freq) {
	return 21.4 * log10(2.71828) * log(4.37*1e-3*_freq + 1.0);
}

void feat_base::rfft(FrmDat *_signal)
{
	int n1, n2, i, j, k, l;
	BaseFloat xt, yt, c, s;
	BaseFloat e, a;
	int m = 0;
	int n = 0;
	BaseFloat *tmp = new BaseFloat[_signal->frmSize];
	while(pow(2.0,m)<_signal->frmSize)
		 m++;
	n = pow(2.0,m);
	for (int i=0; i<_signal->frmSize; i++)
		tmp[i] = _signal->_real[i];

	if (n!=_signal->frmSize) {
		delete []_signal->_real;
		_signal->_real = new BaseFloat[n];
		for (int i=0; i<_signal->frmSize; i++)
			_signal->_real[i] = tmp[i];
		for (int i=_signal->frmSize; i<n; i++)
			_signal->_real[i] = 0.0;
	}
	if (_signal->_img==0)
		_signal->_img = new BaseFloat[n];
	if (_signal->_power==0)
		_signal->_power = new BaseFloat[n];
	for (int i=0; i<n; i++) {
		_signal->_img[i] = 0.0;
		_signal->_power[i] = 0.0;
	}
	_signal->frmSize = n;

	/* Loop through all m stages */
	n2 = n;
	for (k=0; k<m; k++) {
		n1 = n2;
		n2 = n2 / 2;
		e = 2 * M_PI / n1;
		for (j=0; j<n2; j++) {
			/* Compute Twiddle factors */
			a = j * e;
			c = (BaseFloat)cos(a);
			s = (BaseFloat)sin(a);

			/* Do the butterflies */
			for (i=j; i<n; i+=n1) {
				l = i + n2;
				xt = _signal->_real[i] - _signal->_real[l];
				_signal->_real[i] = _signal->_real[i] + _signal->_real[l];
				yt = _signal->_img[i] - _signal->_img[l];
				_signal->_img[i] = _signal->_img[i] + _signal->_img[l];
				_signal->_real[l] = c*xt + s*yt;
				_signal->_img[l] = c*yt - s*xt;
			}
		}
	}
	/* Bit reversal : descrambling */
	j = 0;
	for (int i=0; i<n-1; i++) {
		if (i<j) {
			xt = _signal->_real[j];
			_signal->_real[j] = _signal->_real[i];
			_signal->_real[i] = xt;
			xt = _signal->_img[j];
			_signal->_img[j] = _signal->_img[i];
			_signal->_img[i] = xt;
		}
		k = n / 2;
		while (k<=j) {
			j -= k;
			k /= 2;
		}
		j += k;
	}
	delete []tmp;
}
}
