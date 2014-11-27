#ifndef _LABEL_H
#define _LABEL_H

#include "../base/kaldi-common.h"
#include "../matrix/jun-math.h"

using namespace kaldi;
using namespace std;

/* Reading MLF file and search certain transcription name */
namespace jun {

#define Multiplier -1664117991L
#define Maxstrlen 1024
#define MaxHashBuckets	700
#define INF	FLT_MAX
#define REAL_EPSILON	FLT_EPSILON
#define LOG_ZERO	-INF

/* Structure of word */
struct Word {
	string WordName;
	int16 Id;
	Word *Prev;
	Word *Next;
	explicit Word() : Id(0), Prev(0), Next(0)
	{}
};

/* Structure of transcription name */
struct Trans_label {
	Word *wrd_link;
	string trans_fn;
	int16 wrd_total;
	Trans_label *Prev;
	Trans_label *Next;
	explicit Trans_label() : wrd_link(0), wrd_total(0),
						Prev(0), Next(0)
	{}
};

/* Hash Table */
struct TransHashTable {
	Trans_label *trans_link;
	int num_transFn;
	explicit TransHashTable() : trans_link(0), num_transFn(0)
	{}
};

class Transcription {
private:
	int	Hash(string s);
	TransHashTable *HashLink;
public:
	Transcription() {
		HashLink = new TransHashTable[MaxHashBuckets];
	}
	virtual ~Transcription() {}
	bool Read(string Mlf_fn);
	Trans_label* Search(string Trans_fn);
	void Output(string Trans_fn);
	void Delete();
};

// Structure of HTK header
struct HTKhdr {
    int32 nSamples;         // number of samples
    int32 sampPeriod;       // samples period
    int16 sampSize;         // sample point size (sizeof(short)*nSamples)
    int16 parmKind;         // value of feature parameter kind
    std::string parmName;   // Name of feature parameter
};

void ReadHTKHeader(std::ifstream &in, HTKhdr *header, bool Swap);
bool ParseHTKParmKind(int16 parmKind, std::string &outstr);

}

#endif	// _LABEL_H
