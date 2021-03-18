#! /bin/bash
gcc multigrid_par.c -fopenmp -o multigrid_par.out
echo -n > ./multigrid_par.dat
iterations=$(seq 5)
numWorkers=$(seq 4) 

for w in $numWorkers; do
  size="size: 12"
  for j in $iterations; do
    echo $size "numWorkers: " $w $(./multigrid_par.out 12 12500000 $w) >> ./temp$w.dat
  done

  size="size: 24"
  for n in $iterations; do
  echo $size "numWorkers: " $w $(./multigrid_par.out 24 3000000 $w) >> ./temp$w.dat
  done

  median100=$(sort ./temp$w.dat -g -k "2,11" | awk 'NR==3')
  echo $median100 >> ./multigrid_par.dat
  median200=$(sort ./temp$w.dat -g -k "2,11" | awk 'NR==8')
  echo $median200 >> ./multigrid_par.dat
  rm ./temp$w.dat
done