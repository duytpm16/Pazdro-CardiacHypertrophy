### Parallel permutation scans. Looping from 1-13 because that's how many phenotypes there are in the dataset

for i in {1..13}
do
  echo "#PBS -l nodes=1:ppn=12,mem=16gb
  module load R/3.5.1

  Rscript pazdro_sex_effect_permutation_HPC.R $i 10000 12" >> pazdro_sexeffect_perm_${i}.sh
  qsub pazdro_sexeffect_perm_${i}.sh
done
