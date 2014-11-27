#ifndef JUN_WAVE_H_
#define JUN_WAVE_H_	1

#include <cstdio>
#include <string.h>

using namespace std;

namespace jun {
// ��Ƶ�ṹ�嶨��
struct Audio {
	short* data;		// ��������
	int totLen;			// �ܳ��� ����wav�ļ���С��ȥID��Size��ռ�õ��ֽ���
	int BlockLen;		// �鳤��
	short WavKind;		// ���뷽ʽ
	short chans;		// ������Ŀ
	int nSmpRate;		// ����Ƶ��
	int nBytepSec;		// ÿ�������ֽ���
	short nBytepSmp;	// ÿ��������Ҫ���ֽ���
	short nBitpSmp;		// ÿ��������Ҫ��λ��
	int dataLen;		// ���ݴ洢����
	explicit Audio() : data(0), totLen(0), BlockLen(0), WavKind(1),
		chans(22), nSmpRate(16000), nBytepSec(2), nBytepSmp(4),
		nBitpSmp(16), dataLen(0)
	{}
};

//��Ƶ��
class Wave {
	const char* wavfilename;
	void AddWavHeaderInfo(FILE* _fp, Audio& _audio);	// дͷ�ļ�
	bool ReadWavHeaderInfo(FILE* _fp, Audio& _audio);	// ��ͷ�ļ�
public:
	Wave(const char* _wavfile) : wavfilename(_wavfile) {}
	~Wave() {}
	void ReadWavData(Audio &audio, bool ReadWavHead);						// ������
	void Write2Wav(const char* _outFn, Audio& _audio);	// д����
};
}

#endif ///:~
