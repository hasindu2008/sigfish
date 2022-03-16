#!/bin/bash

read_id=$1

slow5tools get --to slow5 test/sequin_reads.blow5 $1 | grep -v '^[#@]' | awk '{print $8}' > $1.txt &&  matlab.exe -nodisplay -nosplash -nodesktop -r "sig=dlmread('$1.txt'); plot(sig); savefig('$1.fig')"

