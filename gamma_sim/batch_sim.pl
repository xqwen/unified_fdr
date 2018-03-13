$count = 1;
for ($s = 0.1;$s<= 0.9; $s += 0.1){
   for(1..100){	
	printf "Rscript sim.R $s $count\n";
	$count++;
   }
}
