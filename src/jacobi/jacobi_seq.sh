#! /bin/bash
gcc jacobi_seq.c -o jac_seq.out
echo -n > ./jacobi_seq.dat
iterations=$(seq 5)
size="size: 100"
for j in $iterations; do
  echo $size $(./jac_seq.out 100 500000) >> ./temp.dat
done

size="size: 200"
for n in $iterations; do
 echo $size $(./jac_seq.out 200 125000) >> ./temp.dat
done

median100=$(sort ./temp.dat -g -k "2,7" | awk 'NR==3')
echo $median100 > ./jacobi_seq.dat
median200=$(sort ./temp.dat -g -k "2,7" | awk 'NR==8')
echo $median200 >> ./jacobi_seq.dat

rm ./temp.dat