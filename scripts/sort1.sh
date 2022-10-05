
input=$1
subset=$2
tosort=$3
#EXAMPLE
# input=wuhan_vs_alpha_infected_child.tsv.gz; subset=pv_pb; tosort=pv_pb_rand
# module load anaconda3/5.3.0
#sh sort.sh Alpha_infected_vs_uninfected_child.tsv.gz *pv*pb pv_pb_rand | less -S 

#source activate csvkit
##NOTE I SHOULD HAVE CALLED THE ENVIRONMENT CSVTK AS JUST USING THAT TOOL
 pv_lines=$(zcat $input  | head -n 1 | csvtk transpose -t | grep $subset | csvtk transpose)
pv_lines1="ID,${pv_lines}"
echo "#${pv_lines}"
# zcat $1 Alpha_infected_vs_uninfected_child.tsv.gz  | cut -f 1,6,9,10,11 | grep -v '^ENSG'  | csvsort -t -c wilcox_pv_cell | less -S


#zcat Alpha_infected_vs_uninfected_child.tsv.gz | awk 'NF==19{print}{}' | csvtk cut -d '	'  -F -f 'ID,*pv*pb*'   | csvtk sort   -k 'wilcox_pv_pb'  | less -S

if [ $tosort ] ; then
 pvline=$(echo $pv_lines  | csvtk transpose  | grep $tosort | csvtk transpose)
 echo "#${pvline}"
zcat $input | grep -v '^ENSG' | awk 'NF==19{print}{}' | csvtk cut -d '	'  -f $pv_lines1   | csvtk sort   -k $pvline 
elif [ $subset ]; then
zcat $input | grep -v '^ENSG' | awk 'NF==19{print}{}' | csvtk cut -d '	'   -f $pv_lines1
fi	
