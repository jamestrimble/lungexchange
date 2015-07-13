mkdir -p temp
echo 'Generating donors and patients...'
python lungexchange/lung_generator.py 100000 > temp/temp.txt

echo 'Showing distributions...'
NPATIENTS=$(head -n 1 temp/temp.txt | cut -d ' ' -f 1) 
NDONORS=$(head -n 1 temp/temp.txt | cut -d ' ' -f 2)
sed -n "2,$((NPATIENTS+1)) p" temp/temp.txt > temp/patients.txt
sed -n "$((NPATIENTS+2)),$((NPATIENTS+NDONORS+1)) p" temp/temp.txt > temp/donors.txt

python scripts/check_probs.py temp/patients.txt temp/donors.txt

rm temp/temp.txt temp/patients.txt temp/donors.txt
 
