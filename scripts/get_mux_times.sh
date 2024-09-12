#!/bin/bash

# run on: pore_scan_data_XXXXXXX.csv

scan_data=$1
split=0 # for splitting data on odd and even channels

# temp data
freq=freq.csv
entries=entries.csv

die() {
    echo "$@" 1>&2
    exit 1
}

# cut columns: channel=1, mux_scan_asessment=27, seconds_since_start_of_run=37
{ awk -F, '{ print $37" "$27" "$1 }' $scan_data | tail -n +2 | sort -nk1,2 > $entries; } || die "failed to collect mux data entries"

if [ $split -eq 1 ]
then
    { grep "^.*[02468]$" $entries | awk '{ print $1" "$2 }' | uniq -c | awk '{ print $2","$3","$1 }' > "even_${freq}"; } || die "failed to collect frequnecy data"
    { grep "^.*[13579]$" $entries | awk '{ print $1" "$2 }' | uniq -c | awk '{ print $2","$3","$1 }' > "odd_${freq}"; } || die "failed to collect frequnecy data"
else
    { awk '{ print $1" "$2 }' $entries | uniq -c | awk '{ print $2","$3","$1 }' > $freq; } || die "failed to collect frequnecy data"
fi

export split=${split}

python3 -c "

import os
import csv

def run(freq_path, out_path):
    # read frquency csv
    table = dict({
        'multiple': {},
        'other': {},
        'reserved_pore': {},
        'saturated': {},
        'single_pore': {},
        'unavailable': {},
        'zero': {},
    })
    timestamps = set()
    
    with open(freq_path, 'r') as file:
        reader = csv.reader(file, delimiter=',')
        for row in reader:
            time = int(row[0])
            eval = row[1]
            freq = int(row[2])
            
            if time not in timestamps:
                timestamps.add(time)
                table['multiple'][time] = 0
                table['other'][time] = 0
                table['reserved_pore'][time] = 0
                table['saturated'][time] = 0
                table['single_pore'][time] = 0
                table['unavailable'][time] = 0
                table['zero'][time] = 0
                
            table[eval][time] = freq
    
    # print frequency table
    with open(out_path, 'w') as file:
        file.write('time(sec),multiple,other,reserved_pore,saturated,single_pore,unavailable,zero\n')
        
        sorted_ts = sorted(timestamps);
        for time in sorted_ts:
            single_pore = float(table['single_pore'][time])
            file.write(
                str(time)
                + ',' + str(table['multiple'][time])
                + ',' + str(table['other'][time])
                + ',' + str(table['reserved_pore'][time])
                + ',' + str(table['saturated'][time])
                + ',' + str(table['single_pore'][time])
                + ',' + str(table['unavailable'][time])
                + ',' + str(table['zero'][time])
                + '\n'
            )

split = int(os.getenv('split'))

if split == 1:
    run('odd_freq.csv', 'odd_freq_table.csv')
    run('even_freq.csv', 'even_freq_table.csv')
else:
    run('freq.csv', 'freq_table.csv')
"

# cleanup
rm $entries || die "could not remove entries file"
if [ $split -eq 1 ]
then
    rm "even_$freq" || die "could not remove frequency file"
    rm "odd_$freq" || die "could not remove frequency file"
else
    rm "$freq" || die "could not remove frequency file"
fi 


echo "all done!"