
OBJFILES = jun-featbase.o jun-featGfcc.o jun-wave.o jun-mainGfcc.o \
					 ../base/kaldi-base.a ../util/kaldi-util.a
					 
Gfcc_extract : ${OBJFILES}
	g++ -O3 -o Gfcc_extract ${OBJFILES}	

Wave.o : Wave.cc
feat-base.o : jun-featbase.cc
Gfcc_feat.o : jun-featGfcc.cc
Gfcc_main.o : jun-mainGfcc.cc

clean:
	-rm -f *.o *.a $(TESTFILES) tmp* *.tmp

depend:
	-$(CXX) -M $(CXXFLAGS) *.cc > .depend.mk  

# removing automatic making of "depend" as it's quite slow.
#.depend.mk: depend

-include .depend.mk