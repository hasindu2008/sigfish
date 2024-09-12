#!/bin/bash

out="fish_vs_nofish.csv"
fish="odd_freq_table.csv"
nofish="even_freq_table.csv"

die() {
    echo "$@" 1>&2
    exit 1
}

if [ ! -f "$fish" ]
then
    echo "File $fish does not exist"
fi

if [ ! -f "$nofish" ]
then
    echo "File $nofish does not exist"
fi

# get the time  respective columns of each file
echo "time(sec),with_sigfish,without_sigfish" > $out || die "script failed"
{ join -t"," -j 1 -o 1.1,1.6,2.6 $fish $nofish | tail -n +2 >> $out; } || die "script failed"

echo "all done!"