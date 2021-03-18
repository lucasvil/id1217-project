#! /bin/bash
gcc multigrid_seq.c -o multigrid_seq.out
echo -n > ./multigrid_seq.dat
iterations=$(seq 5)
size="size: 12"
for j in $iterations; do
  echo $size $(./multigrid_seq.out 12 12500000) >> ./temp.dat
done

size="size: 24"
for n in $iterations; do
 echo $size $(./multigrid_seq.out 24 3000000) >> ./temp.dat
done

median100=$(sort ./temp.dat -g -k "2,7" | awk 'NR==3')
echo $median100 > ./multigrid_seq.dat
median200=$(sort ./temp.dat -g -k "2,7" | awk 'NR==8')
echo $median200 >> ./multigrid_seq.dat

rm ./temp.dat