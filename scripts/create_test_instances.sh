let END=20
for size in 15 30; do
  for i in $(seq 1 $END); do
    python lungexchange/lung_generator.py $size > test_instances/$size-$i.txt
  done
done
