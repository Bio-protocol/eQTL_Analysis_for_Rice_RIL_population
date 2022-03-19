##command 1 
QTLtools trans --vcf eQTL_genotype.vcf.gz --bed flag_leaf_eTrait.bed.gz --nominal --threshold 0.05 --out flag_leaf005.trans.nominal.hits.txt.gz > running.log

##command 2
for i in {1..100};do
       QTLtools trans --vcf eQTL_genotype.vcf.gz --bed flag_leaf_eTrait.bed.gz --threshold 0.05 --permute --out flag_leaf005_trans_perm_${i} --seed ${i} > running.log
done

##command 3
zcat flag_leaf005_trans_perm_*.hits.txt.gz | gzip -c > flag_leaf005_permutations_all.txt.gz
Rscript runFDR_ftrans.R flag_leaf005.trans.nominal.hits.txt.gz flag_leaf005_permutations_all.txt.gz flag_leaf_trans_005_permutations_all.txt

