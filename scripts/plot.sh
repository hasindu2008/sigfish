#!/bin/bash

set -x
set -e

die () {
    echo >&2 "$@"
    exit 1
}

if [ $# -ne 3 ]; then
    echo "Usage: $0 <read_id> <file.blow5> <file.fastq>"
    exit 1
fi

#scripts/plot.sh 00213403-4297-4f03-8412-3cc8b9cb845a test/sequin_reads.blow5 test/sequin_reads.fastq

blow5=${2}
fastq=${3}
read_id=${1}

rm -f ${read_id}.txt ${read_id}.events.txt ${read_id}.adaptor.txt ${read_id}.fasta ${read_id}.abea.txt ${read_id}.eventalign.txt ${read_id}.events.fig

slow5tools get --to slow5 ${blow5} ${read_id} | grep -v '^[#@]' | awk '{print $8}' > ${read_id}.txt
./sigfish event ${blow5} ${read_id} -n | awk '{print $3"\t"$4"\t"$5}' > ${read_id}.events.txt
./sigfish seg ${blow5} ${read_id} -n | cut -f 3,4,5,6 > ${read_id}.adaptor.txt

samtools faidx ${fastq}
samtools faidx ${fastq} ${read_id} > ${read_id}.fasta
f5c resquiggle --rna ${read_id}.fasta ${blow5} | tail -n +2 | awk '{print $5"\n"$6}' | datamash min 1 max 1 > ${read_id}.abea.txt

#zgrep ${read_id}  test/sequin_reads.abea.tsv.gz | awk '{print $5"\n"$6}' | datamash min 1 max 1 > ${read_id}.abea.txt
#zgrep ${read_id}  test/sequin_reads.eventalign.tsv.gz | awk '{print $NF"\n"$(NF-1)}' | datamash min 1 max 1 > ${read_id}.eventalign.txt
echo "0 0" > ${read_id}.eventalign.txt

matlab.exe -nodisplay -nosplash -nodesktop -minimize -r "
a=dlmread('${read_id}.txt'); b=dlmread('${read_id}.events.txt');
startidx=b(:,1)+1; endidx=b(:,2);
avg=zeros(length(a),1);
for j=1:length(startidx)
    avg(startidx(j):endidx(j))=mean(a(startidx(j):endidx(j)));
end
x=dlmread('${read_id}.adaptor.txt');
x1=x(:,[1:2]);
x2=x(:,[3:4]);
y=[1200,1200];
x3=dlmread('${read_id}.abea.txt');
%x4=dlmread('${read_id}.eventalign.txt');

f=figure;
plot(a); hold on; plot(avg); xlabel('sample index'), ylabel('raw signal value'); stem(x1,y); stem(x2,y); stem(x3,y); legend('raw signal','events','jnn-adaptor','jnn-polya', 'abea', 'eventalign');
savefig(f,'${read_id}.fig');
saveas(f,'${read_id}.png');
%close all;
%c=dlmread('${read_id}.adaptor.txt');
%figure; plot(b(:,3)); savefig('${read_id}.events.fig');
%exit
"


