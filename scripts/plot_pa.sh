#!/bin/bash

#scripts/plot_pa.sh 00213403-4297-4f03-8412-3cc8b9cb845a test/sequin_reads.blow5
# 00213403-4297-4f03-8412-3cc8b9cb845a
# 00acf97d-7374-4587-ae97-a6605281cbbf
# 01046039-3547-46be-bfec-9884cd9ecb7d
# 0139217f-99c5-4773-954f-7441f70b65de
# 017eed00-82de-4a63-997a-b4aeed8fc741
# 01906a38-cf7f-442c-a837-c5097c16cbf1
# 01a585aa-b7d5-4dcc-92d8-a556a4a41df8
# 020efcb0-c412-4f37-a993-7ab3e63248df
# 03239181-a3ce-41ba-869d-98d83b35c177

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <read_id> <file.blow5>"
    exit 1
fi

FILE=${2}
read_id=${1}

#FILE=test/sequin_reads.blow5
#FILE="/mnt/c/Users/hasin/Dropbox (Garvan)/ongoing/sigfish/collapse_ref/250_events_250_ref_analysis_deep/extracted_belong_1_bad_map.blow5"

rm -f ${read_id}.pa.txt

./sigfish pa "$FILE" ${read_id} |  cut -f 3 | tail -n +2 | tr ',' '\t' > ${read_id}.pa.txt

matlab.exe -nodisplay -nosplash -nodesktop -r "
a=dlmread('${read_id}.pa.txt');
plot(a); hold on; xlabel('sample index'), ylabel('raw signal value'); legend('raw signal'); title('$read_id');
"

