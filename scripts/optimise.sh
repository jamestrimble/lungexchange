for f in instances/??-*.txt; do
  python lungexchange/optimise.py $f 5 True
done
for f in instances/??-*.txt; do
  python lungexchange/optimise.py $f 5 False
done
