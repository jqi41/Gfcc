#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "../feat/jun-wave.h"

namespace jun {

const int WAVE_FORMAT_PCM	=	1;
const int WAVE_FORMAT_MULAW	=	0x0006;
const int WAVE_FORMAT_ALAW	=	0x0007;

/*
// RIFF WAVE Chunk
0x52�� 0x49�� 0x46�� 0x46�� // "RIFF"
0x30�� 0x00�� 0x00�� 0x00�� // �ܳ��� ����wav�ļ���С��ȥID��Size��ռ�õ��ֽ���
0x57�� 0x41�� 0x56�� 0x45�� // "WAVE"
// Format Chunk
0x66�� 0x6D�� 0x74�� 0x20�� // "fmt "
0x10�� 0x00�� 0x00�� 0x00�� // �鳤��
0x01�� 0x00�� // ���뷽ʽ wFormatTag
0x01�� 0x00�� // ������Ŀ wChannels
0x80�� 0x3E�� 0x00�� 0x00�� // ����Ƶ�� dwSamplesPerSec
0x00�� 0x7D�� 0x00�� 0x00�� // ÿ�������ֽ��� dwAvgBytesPerSec
0x02�� 0x00�� // ÿ��������Ҫ���ֽ��� wBlockAlign
0x10�� 0x00�� // ÿ��������Ҫ��λ�� wBitsPerSample
// Data Chunk
0x64�� 0x61�� 0x74�� 0x61�� // "data"
0x00�� 0x00�� 0x00�� 0x00�� // �鳤��
*/
bool Wave::ReadWavHeaderInfo(FILE* pcm_fp, Audio &audio) {
	char magic[4];

	fseek(pcm_fp, 0, SEEK_SET);
	fread(magic, 4, 1, pcm_fp);
	if (strncmp("RIFF", magic, 4)) {
		cout << "Input file is not in RIFF format." << endl;
		return false;
	}
	fread(&audio.totLen, 4, 1, pcm_fp);
	fread(magic, 4, 1, pcm_fp);
	if (strncmp("WAVE", magic, 4)) {
		cout << "Input file is not in WAVE format" << endl;
		return false;
	}

	/* Looking for "fmt" before the end of file */
	while(1) {
		if (feof(pcm_fp)) {
			cout << "No data portion in WAVE file" << endl;
			return false;
		}
		fread(magic, 4, 1, pcm_fp);
		if (strncmp("data", magic, 4)==0) {
			fread(&audio.dataLen, 4, 1, pcm_fp);
			break;
		}
		else {
			fread(&audio.BlockLen, 4, 1, pcm_fp);
		}
		if (strncmp("fmt ", magic, 4)==0)
		{
			fread(&audio.WavKind, 2, 1, pcm_fp);
			if (audio.WavKind != WAVE_FORMAT_PCM)
			{
				cout << "Only standard PCM is supported" << endl;
				return false;
			}
			fread(&audio.chans, 2, 1, pcm_fp);
			if (audio.chans!=1)
			{
				cout << "No single channel" << endl;
				return false;
			}
			fread(&audio.nSmpRate, 4, 1, pcm_fp);
			if (audio.nSmpRate!=16000 && audio.nSmpRate!=8000) {
				cout << "Sample rate" << audio.nSmpRate << " is not supported" << endl;
				return false;
			}
			fread(&audio.nBytepSec, 4, 1, pcm_fp);
			fread(&audio.nBytepSmp, 2, 1, pcm_fp);
			fread(&audio.nBitpSmp, 2, 1, pcm_fp);
			if (audio.nBitpSmp!=16 && audio.nBitpSmp!=8)
			{
				cout << "Only 8/16 bits audio supported" << endl;
				return false;
			}
			if (audio.WavKind==WAVE_FORMAT_MULAW || audio.WavKind==WAVE_FORMAT_ALAW) {
				if (audio.nBitpSmp != 8) {
					cout << "Only 8-bit mu-law/a-law is supported" << endl;
					return false;
				}
				if (audio.nSmpRate != 8000) {
					cout << "Only 8k mu-law/a-law is supported" << endl;
					return false;
				}
			}
		}
	}
	return true;
}

void Wave::ReadWavData(Audio &audio, bool ReadWavHead)
{
	FILE *pcm_fp = 0;
	if ((pcm_fp=(fopen(wavfilename, "rb")))==NULL)	
		exit(-1);
	if (ReadWavHead)
		ReadWavHeaderInfo(pcm_fp, audio);
	int nSize;
	fseek(pcm_fp, 0, SEEK_END);
	
	if (ReadWavHead) 
		nSize = ftell(pcm_fp) - 44;
	else {
		nSize = ftell(pcm_fp);
		audio.dataLen = nSize;
	}
	nSize /= sizeof(short);
	audio.data = new short[nSize];
	if (audio.data==NULL)
		audio.data = new short[nSize];
	memset(audio.data, 0, nSize);
	if (ReadWavHead)
		fseek(pcm_fp, 44, SEEK_SET);
	else
		fseek(pcm_fp, 0, SEEK_SET);
	fread(audio.data, sizeof(short), nSize, pcm_fp);
	fclose(pcm_fp);
	return;
}

void Wave::AddWavHeaderInfo(FILE *wp, Audio& audio)
{
	fseek(wp, 0, SEEK_SET);
	fwrite("RIFF", 4, 1, wp);	fwrite(&audio.totLen, 4, 1, wp);
	fwrite("WAVE", 4, 1, wp);
	fwrite("fmt ", 4, 1, wp);	fwrite(&audio.BlockLen, 4, 1, wp);
	fwrite(&audio.WavKind, 2, 1, wp);
	fwrite(&audio.chans, 2, 1, wp);
	fwrite(&audio.nSmpRate, 4, 1, wp);
	fwrite(&audio.nBytepSec, 4, 1, wp);
	fwrite(&audio.nBytepSmp, 2, 1, wp);
	fwrite(&audio.nBitpSmp, 2, 1, wp);
	fwrite("data", 4, 1, wp);
	fwrite(&audio.dataLen, 4, 1, wp);
}

void Wave::Write2Wav(const char* outFn, Audio& audio)
{
	FILE *out_fp = 0;
	if ((out_fp=fopen(outFn, "wb"))==NULL)
		exit(-1);
	AddWavHeaderInfo(out_fp, audio);
	fwrite(audio.data, sizeof(short), audio.dataLen / sizeof(short), out_fp);

	fclose(out_fp);
	return;
}
}
