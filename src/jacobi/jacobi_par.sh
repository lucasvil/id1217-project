#! /bin/bash
gcc jacobi_par.c -fopenmp -o jac_par.out
echo -n > ./jacobi_par.dat
iterations=$(seq 5)
numWorkers=$(seq 4) 

for w in $numWorkers; do
  size="size: 100"
  for j in $iterations; do
    echo $size "numWorkers: " $w $(./jac_par.out 100 500000 $w) >> ./temp$w.dat
  done

  size="size: 200"
  for n in $iterations; do
  echo $size "numWorkers: " $w $(./jac_par.out 200 125000 $w) >> ./temp$w.dat
  done

  median100=$(sort ./temp$w.dat -g -k "2,9" | awk 'NR==3')
  echo $median100 >> ./jacobi_par.dat
  median200=$(sort ./temp$w.dat -g -k "2,9" | awk 'NR==8')
  echo $median200 >> ./jacobi_par.dat
  rm ./temp$w.dat
done