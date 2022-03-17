#!/bin/bash

read_id=$1

slow5tools get --to slow5 test/sequin_reads.blow5 $1 | grep -v '^[#@]' | awk '{print $8}' > $1.txt 
./sigfish events --rna test/sequin_reads.blow5 $1 | tail -n +2 | awk '{print $3"\t"$4"\t"$5}' > $1.events.txt 
matlab.exe -nodisplay -nosplash -nodesktop -r "
a=dlmread('$1.txt'); b=dlmread('$1.events.txt');
startidx=b(:,1)+1; endidx=b(:,2)+startidx-1; 
avg=zeros(length(a),1);
for j=1:length(startidx)
    avg(startidx(j):endidx(j))=mean(a(startidx(j):endidx(j)));
end
plot(a); hold on; plot(avg); savefig('$1.fig') 
"

