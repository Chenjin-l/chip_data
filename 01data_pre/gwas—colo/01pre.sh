
awk 'BEGIN{OFS="\t"}
NR==1{next}
$7 < 5e-8{
    chr=$1;
    pos=$3;
    snp=$2;
    print chr, pos-1, pos, snp, $7
}' gwas_eye_muscle_step1.step2.assoc > eye_muscle.sig.snp.bed

bedtools intersect \
-a eye_muscle.sig.snp.bed \
-b /home/Mzhou/zPig_data/04Pig_Gtex/01molQTL_Data/PigGTEx_v0.permutations_eQTL/eqtl_bed/*.eqtl.bed \
-wa -wb > eye_muscle.eqtl.overlap.txt
