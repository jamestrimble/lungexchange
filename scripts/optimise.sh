for f in instances/??-1??.txt; do
  python lungexchange/optimise.py $f 5 True
done
for f in instances/??-1??.txt; do
  python lungexchange/optimise.py $f 5 False
done
