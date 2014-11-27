#!/bin/bash

echo "make kaldi-base.a"
cd base 
make 
cd ..

echo "make kaldi-util.a"
cd util
make
cd ..

echo "Generate gfcc_extract"
cd feat
make
cd ..
cd ..

mv ./src/feat/Gfcc_extract ./bin/Gfcc_extract