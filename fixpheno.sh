#!/bin/sh


while (("$#")); do
	if [ -e $1.pheno.txt ]; then
		echo Fixing $1.pheno.txt;
            	head -1 $1.pheno.txt > tmp.pheno.txt;
            	tail -n+2 $1.pheno.txt | sed -r 's/ /\./g;s/\//\./g;s/(NULL|NA|#NUM!|-Inf|Inf)/-9/g' >> tmp.pheno.txt;
            	mv -f tmp.pheno.txt $1.pheno.txt;
        fi;
        if [ -e $1.covar.txt ]; then
            	echo Fixing $1.covar.txt;
            	sed -i 's/ /\./g;s/\//\./g;s/(NULL|NA|#NUM!|-Inf|Inf)/-9/g' $1.covar.txt;
        fi;
        shift;
done
