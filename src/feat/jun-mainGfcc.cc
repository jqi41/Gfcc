#include "../feat/jun-featGfcc.h"
#include "../util/parse-options.h"
#include <time.h>

using namespace std;
using namespace kaldi;
using namespace jun;

/* Space is used to split one line into speech and feature search paths */
bool ParseUnsSym(std::string &buf, int16 &init, int16 &st, int16 &ed) {
	unsigned int i=0, j=0, m=0;
	bool space = false;
	init = st = ed = 0;
	while(i<buf.size()) {
		if(isspace(buf.at(i))==0) {
			init = i;		
			break;
		}
		i++;
	}
	if (buf.at(init)=='#')
		return false;

	i = 0;
	while (i<buf.size()) {
		if (isspace(buf.at(i))!=0 && space==false) {
			ed = i;
			space = true;
		}
		else if (isspace(buf.at(i))==0 && space==true) {
			st = i;
			space = false;
		}
		i++;
	}
	return true;
}

int main(int argc, char *argv[]) {
    Audio audio;
	clock_t start, finish;
	BaseFloat duration;
	const char *usage = 
		"Create GFCC feature files.\n"
		"Single Mode Usage: Gfcc_extract [options...] <wav-rspecifier> <feats-wspecifier>\n"
		"Batch Mode Usage: Gfcc_extract [options...] <scp_rspecifier>\n";
	ParseOptions po(usage);
	Gfcc_Options gfcc_opt;

	gfcc_opt.Register(&po);
	feat_gfcc gfcc;
	po.Read(argc, argv);
	int16 numArg = po.NumArgs();
	if (numArg != 2 && numArg != 1) {
		po.PrintUsage();
		exit(1);
	}
	if (numArg == 2) {
		start = clock();
		gfcc_opt.WavFn.append(po.GetArg(1));
		gfcc_opt.FeatFn.append(po.GetArg(2));
		gfcc.Gfcc_Gen(gfcc_opt);
		finish = clock();
		duration = (BaseFloat)(finish - start) / CLOCKS_PER_SEC;
		std::cout << duration << std::endl;
	}
	else if (numArg == 1) {
		std::fstream file((po.GetArg(1)).c_str(), std::ios::in);
		if (file.fail()) {
			std::cout << "Cannot open coding search file !" << std::endl;
			exit(1);
		}
		int16 st=0, ed=0, init=0;
		string line, tmp;
		char buf[MAXSTRLEN] = "";
		while (file.good()) {
			line.clear();
			file.getline(buf, MAXSTRLEN);
			line.append(buf);
			if (line.empty()) 
				continue;
			if (!ParseUnsSym(line, init, st, ed))
				continue;
			tmp.clear();
			for (size_t k=st; k<line.size(); k++)
				tmp += buf[k];
			gfcc_opt.FeatFn.append(tmp.c_str());
			tmp.clear();
			for (int16 k=init; k<ed; k++)
				tmp += buf[k];
			gfcc_opt.WavFn.append(tmp.c_str());
			tmp.clear();
			gfcc.Gfcc_Gen(gfcc_opt);
			gfcc_opt.FeatFn.clear();
			gfcc_opt.WavFn.clear();
		}
		line.clear();
	}
	return 0;
}