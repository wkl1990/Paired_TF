cd 01.rawdata
rename 's/_S[0-9]{1,3}[_L0-9]*_R/_R/' ${1}*.fastq.gz
rename 's/_001.fastq/.fq/' ${1}*.fastq.gz
ls -alt | sed 's/ /\t/g'| sed 's/\t\t/\t/g'|cut -f10| sed 's/_/\t/g'| sort| uniq| cut -f1| sort| uniq | sed ":a;N;s/\n/\t/g;ta"
