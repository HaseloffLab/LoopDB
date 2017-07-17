while read line; do python ../to_fasta.py $line >> fasta.txt; done < list
