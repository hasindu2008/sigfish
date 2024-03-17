# Test Workflows

## Checking accuracy with [Squigulator](https://github.com/hasindu2008/squigulator)

### R9

```sh
# r9 simulate
squigulator -x rna-r9-prom test/rnasequin_sequences_2.4.fa -o rna.blow5 -c rna.paf --prefix=yes --paf-ref

# static dtw
./sigfish dtw test/rnasequin_sequences_2.4.fa rna.blow5 --rna -p -1 -q 500  > dtw.paf
./sigfish eval rna.paf dtw.paf --tid-only

# real-time dtw
./sigfish real test/rnasequin_sequences_2.4.fa rna.blow5 --rna -p -1 -q 500  -t16 > dtw.paf
./sigfish eval rna.paf dtw.paf --tid-only
```

### RNA004

```sh
# rna004 simulate
squigulator -x rna004-prom test/rnasequin_sequences_2.4.fa -o rna.blow5 -c rna.paf --prefix=yes --paf-ref

# static dtw
./sigfish dtw test/rnasequin_sequences_2.4.fa rna.blow5 --rna -p -1 -q 500  > dtw.paf
./sigfish eval rna.paf dtw.paf --tid-only

# real-time dtw
./sigfish real test/rnasequin_sequences_2.4.fa rna.blow5 --rna -p -1 -q 500  -t16 > dtw.paf
./sigfish eval rna.paf dtw.paf --tid-only
```
