#ifndef JUN_WAVE_H_
#define JUN_WAVE_H_	1

#include <cstdio>
#include <string.h>

using namespace std;

namespace jun {
// 音频结构体定义
struct Audio {
	short* data;		// 语音数据
	int totLen;			// 总长度 整个wav文件大小减去ID和Size所占用的字节数
	int BlockLen;		// 块长度
	short WavKind;		// 编码方式
	short chans;		// 声道数目
	int nSmpRate;		// 采样频率
	int nBytepSec;		// 每秒所需字节数
	short nBytepSmp;	// 每个样本需要的字节数
	short nBitpSmp;		// 每个样本需要的位数
	int dataLen;		// 数据存储长度
	explicit Audio() : data(0), totLen(0), BlockLen(0), WavKind(1),
		chans(22), nSmpRate(16000), nBytepSec(2), nBytepSmp(4),
		nBitpSmp(16), dataLen(0)
	{}
};

//音频类
class Wave {
	const char* wavfilename;
	void AddWavHeaderInfo(FILE* _fp, Audio& _audio);	// 写头文件
	bool ReadWavHeaderInfo(FILE* _fp, Audio& _audio);	// 读头文件
public:
	Wave(const char* _wavfile) : wavfilename(_wavfile) {}
	~Wave() {}
	void ReadWavData(Audio &audio, bool ReadWavHead);						// 读数据
	void Write2Wav(const char* _outFn, Audio& _audio);	// 写数据
};
}

#endif ///:~
