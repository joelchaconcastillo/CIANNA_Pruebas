for nday in 20 40 60;
do
  for px in 0.2 0.6;
 do
# for ITV in 0.2;
for ITV in 0.2 0.4; 
  do
    echo "ITV  ".$ITV." Px "${px}
    rm $nday
   for t in 1440 2880;
   #for t in 2880;
   do
	   echo $t
     for i in ../History/*days_${nday}_min_${t}_*_Px_${px}_${ITV};
     do
       tail -10 $i| awk '{print $NF}'| sort -n | head -n 1 >> $nday
     done
     echo "Dias.. $nday"
   Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
   	-e 'sprintf(summary(d), fmt="%.3f")' < $nday
   	#-e 'print(summary(d), digits=20)' < $nday
   done
  done
  done
done
