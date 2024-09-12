#!/bin/bash

log_dir="/mnt/data/dev/logs/minifridge"
log="${log_dir}/${1}.log"
out="${1}.timings"
set -e

die () {
    echo >&2 "$@"
    exit 1
}

total=0
sum=0

shopt -s lastpipe
grep -F 'queries processed in' $log | cut -d' ' -f3,7 | while read -r line;
do
    time_taken=$(echo $line | cut -d' ' -f2)
    sum=$(echo "scale=5; $sum + ($time_taken)" | bc)
    total=$((total+1))
done


echo $total
echo $sum

avg=$(echo "scale=5; $sum / $total" | bc)
echo "average time taken per query: ${avg}"

echo "all done!"
