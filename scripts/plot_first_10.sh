#!/bin/bash

test/eval_rna.sh
grep "rf:Z:tp" err.paf  | awk '{print $1}' > correct.list
grep "rf:Z:fp" err.paf  | awk '{print $1}' > wrong.list

rm *.txt *.fig *.png   -f
head -10 correct.list | while read p; do echo $p; scripts/plot.sh $p; sleep 2; done
rm -f correct/*
mv *.txt *.fig *.png correct/

rm -f wrong/*
rm *.txt *.fig *.png -f
head -10 wrong.list | while read p; do echo $p; scripts/plot.sh $p; sleep 5;  done
mv *.txt *.fig *.png wrong/

