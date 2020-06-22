#! /bin/bash

set -e
# cyfate version
VERSION="0.1"

# default arguments
STHREADS=1
FTHREADS=1
NCYCLE=30
SAMPLE=sample
GENE=gene

# functions for version and help
print_version(){
  echo "cyFate version $VERSION using" `fate.pl 2>&1 1>/dev/null | awk 'NR==1 { print $0 }'`
}

print_help(){
  cat <<EOS
Usage:
cyFate [options] <genome.fna> <query.fa> <protein.faa> <keyword>

cyFate options:
--sthreads INT     Number of threads for search
--fthreads INT     Number of threads for filter
--gene STR         Gene name
--sample STR       Sample name

Extra arguments (-g, -v, etc) will be passed through

Other:
-v|--version         show version info
-h|--help            program usage

For more information, see fate.pl
EOS
}

# parse arguments
while [[ $# -gt 0 ]]; do
  key=${1}
  case ${key} in
    -v|--version)
      print_version
      exit 0
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    --cycle)
      NCYCLE=${2}
      shift
      ;;
    --sthreads)
      STHREADS=${2}
      shift
      ;;
    --fthreads)
      STHREADS=${2}
      shift
      ;;
    --sample)
      SAMPLE="${2}"
      shift
      ;;
    --gene)
      GENE="${2}"
      shift
      ;;
    -*|--*)
      OPTIONS="${OPTIONS} ${1}"
      # OPTIONS="${OPTIONS} ${1} ${2}"
      ;;
    *)
      OPTIONS="${OPTIONS} ${1}"
  esac
  shift
done

SUBJECT=${OPTIONS[0]}
QUERY=${OPTIONS[1]}
PROTEIN=${OPTIONS[2]}
KEYWORD=${OPTIONS[3]}

# echo $subject
echo "cyFate of ${SAMPLE} start"
mkdir -p ${SAMPLE}_cyFate
cd ${SAMPLE}_cyFate

output=${SAMPLE}_${GENE}_i5000_o60_all

for ((ncycle=1; ncycle<${NCYCLE}+1; ncycle++)); do

	echo "${SAMPLE} cycle${ncycle} is started at `date '+%y/%m/%d %H:%M:%S'`"
	mkdir -p cycle_${ncycle}
	cd cycle_${ncycle}

	if [ ${ncycle} = 1 ]; then
		fate.pl search -p ${STHREADS} -h tblastn -g genewise -i 5000 -o 60 -v 1 -x -s ${SUBJECT} ${QUERY} > ${output}.bed || exit
	else
		fate.pl search -p ${STHREADS} -h tblastn -g genewise -i 5000 -o 60 -v 1 -x -s ${SUBJECT} ../${output}_cycle$((ncycle-1))_filter_blue_aa.fasta > ${output}.bed || exit
	fi

	# awk '$9~/0,0,255/{print$0}' ${output}.bed > ${output}_blue.bed
	# bedtools getfasta -fi ${subject} -bed ${output}_blue.bed -s -split -fo ${output}_blue.fasta

	fate.pl filter -p ${FTHREADS} -d ${PROTEIN} -x -h blastx -k "${KEYWORD}" -r 10 ${SUBJECT} ${output}.bed > ${output}_filter.bed || exit
	awk '$9~/0,0,255/{print$0}' ${output}_filter.bed > ${output}_filter_blue.bed
	bedtools getfasta -fi ${SUBJECT} -bed ${output}_filter_blue.bed -s -split -fo ${output}_filter_blue.fasta
		
	translate.py ${output}_filter_blue.fasta > ../${output}_cycle${ncycle}_filter_blue_aa.fasta

	rm -r nucl
	rm -r prot
	echo "${SAMPLE} cycle${ncycle} is finished"
	cd ../
	if [ ${ncycle} == 2 ]; then
		rm ./${output}_cycle$((ncycle-1))_filter_blue_aa.fasta
	elif [[ ${ncycle} != 1 ]]; then
		rm -r ./cycle_$((ncycle-1))
		rm ./${output}_cycle$((ncycle-1))_filter_blue_aa.fasta
	fi
done

cd ../
echo "cyFate of ${SAMPLE} finish"