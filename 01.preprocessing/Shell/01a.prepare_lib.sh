ref=${1}
cat sample.list | while read id || [[ -n ${id} ]]
do
dna=`echo "$id" | cut -f1`
rna=`echo "$id" | cut -f2`
if [[ $dna != "DNA" ]]; then
echo $dna $rna
echo -e "${dna}\t${rna}\t${ref}\t${dna}_${rna}" >> QC_library_list.xls 
fi
done

