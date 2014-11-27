#include "jun-label.h"

namespace jun {
int Transcription::Hash(string s) {
	size_t i;
	unsigned long hashcode = 0;

	hashcode = 01;
	for (i=0; i<s.length(); i++) {
		hashcode = hashcode * Multiplier + s[i];
	}
	return (hashcode % MaxHashBuckets);
}

bool Transcription::Read(string Mlf_fn) {
	if (Mlf_fn.c_str()==0) {
		cout << "Transcription MLF file should be initialized at first."
			<< endl;
		return false;
	}
	ifstream *in = new ifstream;
	in->open(Mlf_fn.c_str(), std::ios::in);

	if (in->fail()) {
		std::cout << "The input stream could not be created !" << std::endl;
		in->close();
		delete in;
		return false;
	}
	in->seekg(0, std::ios::beg);

	char buf[Maxstrlen] = "";
	string labfn;
	string line;
	int16 pos_start=0;
	int16 pos_end=0;
	int Hash_cp;
	TransHashTable *trans_p = 0;

	while (in->good()) {
		line.clear();
		in->getline(buf, Maxstrlen);
		line.append(buf);
		if (line.compare("#!MLF!#")==0) {
			continue;
		}
		else if (line.find(".lab") < line.length()) {
			pos_start = line.find('/');
			pos_end = line.find('.');
			for (int i=pos_start+1; i<pos_end; i++)
				labfn.push_back(buf[i]);

			Hash_cp = Hash(labfn.c_str());
			trans_p = &HashLink[Hash_cp];
			if (trans_p->trans_link==0) {
				trans_p->trans_link = new Trans_label;
				trans_p->trans_link->wrd_link = new Word;
			}
			else {
				while (trans_p->trans_link->Next!=0) {
					trans_p->trans_link = trans_p->trans_link->Next;
				}
				trans_p->trans_link->Next = new Trans_label;
				trans_p->trans_link->Next->Prev = trans_p->trans_link;
				trans_p->trans_link = trans_p->trans_link->Next;
				trans_p->trans_link->wrd_link = new Word;
			}
			trans_p->trans_link->trans_fn.append(labfn.c_str());
			trans_p->num_transFn++;
			labfn.clear();
			continue;
		}
		else if (line.compare(".")==0) {
			trans_p->trans_link->Next = 0;
			trans_p->trans_link->wrd_link->Next=0;
			Hash_cp = 0;
			continue;
		}
		else {
			trans_p->trans_link->wrd_total++;
			trans_p->trans_link->wrd_link->WordName.append(line.c_str());
			trans_p->trans_link->wrd_link->Id = trans_p->trans_link->wrd_total;

			trans_p->trans_link->wrd_link->Next = new Word;
			trans_p->trans_link->wrd_link->Next->Prev = trans_p->trans_link->wrd_link;
			trans_p->trans_link->wrd_link = trans_p->trans_link->wrd_link->Next;
		}
	}
	trans_p = 0;
	in->close();
	delete in;
	return true;
}

Trans_label* Transcription::Search(string trans_fn) {
	int Hash_code=0;
	int idx=0;
	TransHashTable *trans_p=0;

	Hash_code = Hash(trans_fn);
	trans_p = &HashLink[Hash_code];
	while (idx < trans_p->num_transFn) {
		if (trans_fn.compare(trans_p->trans_link->trans_fn.c_str())==0) {
			trans_p->trans_link->Prev = 0;
			trans_p->trans_link->Next = 0;
			return trans_p->trans_link;
		}
		else {
			trans_p->trans_link = trans_p->trans_link->Prev;
			idx++;
		}
	}
	return 0;
}

void Transcription::Output(string Trans_fn) {
	Trans_label *Search_Trans;
	Search_Trans = Search(Trans_fn.c_str());
	if (Search_Trans==0) {
		cout << "The label file " << Trans_fn.c_str() << " cannot be found in the given transcript."
			 << endl;
		return;
	}
	Word p;
	int idx = 0;
	p = *Search_Trans->wrd_link;
	while(idx < Search_Trans->wrd_total) {
		p = *p.Prev;
		idx++;
	}
	idx=0;
	while(idx < Search_Trans->wrd_total) {
		cout << p.WordName.c_str() << " ";
		p = *p.Next;
		idx++;
	}

	cout << endl;
	return;
}

void Transcription::Delete() {
	if (HashLink==0)
		return;
	TransHashTable *tran_p=0;
	for (int i=0; i<MaxHashBuckets; i++) {
		tran_p = &HashLink[i];
		if (tran_p->num_transFn==0)
			continue;
		while (tran_p->trans_link->Next!=0) {
			tran_p->trans_link = tran_p->trans_link->Next;
		}
		while (tran_p->trans_link!=0) {
			while(tran_p->trans_link->wrd_link->Next!=0) {
				tran_p->trans_link->wrd_link = tran_p->trans_link->wrd_link->Next;
			}
			while(tran_p->trans_link->wrd_link!=0) {
				tran_p->trans_link->wrd_link->WordName.clear();
				tran_p->trans_link->wrd_link->Id = 0;
				if (tran_p->trans_link->wrd_link->Prev!=0) {
					tran_p->trans_link->wrd_link = tran_p->trans_link->wrd_link->Prev;
					tran_p->trans_link->wrd_link->Next->Next = 0;
					tran_p->trans_link->wrd_link->Next->Prev = 0;
					delete tran_p->trans_link->wrd_link->Next;
					tran_p->trans_link->wrd_link = 0;
				}
				else {
					tran_p->trans_link->wrd_link->Next = 0;
					tran_p->trans_link->wrd_link->Prev = 0;
					delete tran_p->trans_link->wrd_link;
					tran_p->trans_link->wrd_link = 0;
				}
			}
			if (tran_p->trans_link->Prev!=0) {
				tran_p->trans_link = tran_p->trans_link->Prev;
				tran_p->trans_link->Next->Next = 0;
				tran_p->trans_link->Next->Prev = 0;
				tran_p->trans_link->Next->trans_fn.clear();
				delete tran_p->trans_link->Next;
				tran_p->trans_link->Next = 0;
			}
			else {
				tran_p->trans_link->Next = 0;
				tran_p->trans_link->Prev = 0;
				tran_p->trans_link->trans_fn.clear();
				delete tran_p->trans_link;
				tran_p->trans_link = 0;
			}
		}
	}
	tran_p = 0;	delete tran_p;
	delete [] HashLink;
	HashLink = 0;
	return;
}

bool ParseHTKParmKind(int16 parmKind, std::string &outstr)
{
	int32 tmp = 0;
	tmp = parmKind & 0x003f;
	switch (tmp) {
		case 0:
			outstr.append("WAVEFORM");	break;
		case 1:
			outstr.append("LPC");		break;
		case 2:
			outstr.append("LPREFC");	break;
		case 3:
			outstr.append("LPCEPSTRA");	break;
		case 4:
			outstr.append("LPDELCEP");	break;
		case 5:
			outstr.append("IREFC");		break;
		case 6:
			outstr.append("MFCC");		break;
		case 7:
			outstr.append("FBANK");		break;
		case 8:
			outstr.append("MELSPEC");	break;
		case 9:
			outstr.append("USER");		break;
		case 10:
			outstr.append("DISCRETE");	break;
		case 11:
			outstr.append("PLP");		break;
		default:
			return false;
	}
	tmp = parmKind & 0xffc0;			/* Bit-wise Encoding Qualifier */
	if (tmp & 000100)	outstr.append("_E");
	if (tmp & 000200)	outstr.append("_N");
	if (tmp & 000400)	outstr.append("_D");
	if (tmp & 001000)	outstr.append("_A");
	if (tmp & 002000)	outstr.append("_C");
	if (tmp & 004000)	outstr.append("_Z");
	if (tmp & 010000)	outstr.append("_K");
	if (tmp & 020000)	outstr.append("_0");

	return true;
}

void ReadHTKHeader(std::ifstream &in, HTKhdr *header, bool Swap) {
	int16 sampKind = 0;
	int16 sampSize = 0;
	int32 nSamples = 0;
	int32 sampPeriod = 0;

	in.read((char*)&nSamples, sizeof(nSamples));
	in.read((char*)&sampPeriod, sizeof(sampPeriod));
	in.read((char*)&sampSize, sizeof(sampSize));
	in.read((char*)&sampKind, sizeof(sampKind));

	if ( Swap ) {
		Swap32(&nSamples);
		Swap32(&sampPeriod);
		Swap16(&sampSize);
		Swap16(&sampKind);
	}
	header->nSamples = nSamples;
	header->sampPeriod = sampPeriod;
	header->sampSize = sampSize;
	header->parmKind = sampKind;
	ParseHTKParmKind(header->parmKind, header->parmName);

	return;
}

}
