let END=500
for size in 10 20 50; do
  for i in $(seq 1 $END); do
    python lungexchange/lung_generator.py $size > instances/$size-$i.txt
  done
done
