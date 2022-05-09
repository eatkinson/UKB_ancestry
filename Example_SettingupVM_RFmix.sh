
#reinstall programs on new VM
sudo apt-get update
HOME=/home/eatkinso

PKGS="bzip2 build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libpcre3-dev gfortran openjdk-8-jdk pkg-config"
for p in $PKGS; do
    sudo apt-get install -y $p
done

# Install git.
sudo apt-get install git
cd $HOME

#install Tractor
git clone https://github.com/eatkinson/Tractor.git

#instal BCFtools
wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2
tar -xvf bcftools-1.2.tar.bz2
#ran out of memory. Probably will need to load in a persistent disk to do the actual analysis.
#deleted the huge old reference panel files that were in there

#install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2
tar xvf bcftools-1.2.tar.bz2
cd bcftools-1.2
make
make install
cd htslib-1.2.1
make
make install

#test out bgzip - it lives here: /home/eatkinso/bcftools-1.2/htslib-1.2.1/bgzip
/home/eatkinso/bcftools-1.2/htslib-1.2.1/bgzip temp


#install rfmix
git clone https://github.com/slowkoni/rfmix.git
sudo apt-get install autoconf
sudo apt-get install make
sudo apt install g++

autoreconf --force --install # creates the configure script and all its dependencies
./configure                  # generates the Makefile
make

/home/eatkinso/PGC/rfmix/rfmix --help
#works!


#set up alias paths for ease
alias 'l=ls -lhtr'
alias 'bgzip=/home/eatkinso/bcftools-1.2/htslib-1.2.1/bgzip'
alias 'tabix=/home/eatkinso/bcftools-1.2/htslib-1.2.1/tabix'
alias 'rfmix=/home/eatkinso/PGC/rfmix/rfmix'
alias 'bcftools=/home/eatkinso/PGC/bcftools-1.2/bcftools'


#copy reference files into VM. 

gsutil cp gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schizophrenia/PGC_SCZ_AA_pheno_covs.txt $HOME/PGC/

gsutil cp gs://ukb-diverse-pops/AdmixedAfrEur/Reference-AFR-EUR/HapMapcomb* $HOME/PGC/Afr-Eur-Ref/

gsutil cp gs://ukb-diverse-pops/AdmixedAfrEur/Reference-AFR-EUR/AFR_EUR.1kg_order.indivs.pops.txt $HOME/PGC/Afr-Eur-Ref/

gsutil cp gs://ukb-diverse-pops/AdmixedAfrEur/Reference-AFR-EUR/AFR-EUR-hg19/* $HOME/PGC/Afr-Eur-Ref/
#copied in the recomb map, newly filtered reference panels in gzipped format with tabix index files


##now start processing the VCF files 
zcat vcfheader.txt.gz gpc_mega.cogs_aa.dose.QC.recode.vcf.gz | bgzip > gpc_mega.cogs_aa.dose.QC_1.recode.vcf.gz

# index again and check if bcftools can read it.
tabix -p vcf gpc_mega.cogs_aa.dose.QC_1.recode.vcf.gz
 
/home/eatkinso/PGC/bcftools-1.2/bcftools annotate -x INFO,^FORMAT/GT gpc_mega.cogs_aa.dose.QC_1.recode.vcf.gz > gpc_mega.cogs_aa.dose.QC_2.recode.vcf

bgzip test.vcf
tabix -p vcf test.vcf.gz
/home/eatkinso/PGC/bcftools-1.2/bcftools annotate -x INFO,^FORMAT/GT test.vcf.gz > test2.vcf

#see if works on the first chunk of this file
tabix -p vcf test.vcf.gz

## using the -x FORMAT flag removes all format fields but GT
/home/eatkinso/PGC/bcftools-1.2/bcftools annotate -x FORMAT test.vcf.gz > test.vcf

bgzip test.vcf
tabix -p vcf test.vcf.gz

#also confirm MAF filter works on the test file
bcftools view -q 0.01:minor test.vcf.gz > test1.vcf.gz


....


##then run RFmix
#turn off EM and don't reanalyze reference for speed.
for i in {1..22}; do rfmix -f PGC-SCZ.Bigdeli_merged.aa.MAF1.vcf -r /home/eatkinso/PGC/Afr-Eur-Ref/AFR_EUR.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.maf05.chr$i.recode.vcf.gz --chromosome=${i} -m /home/eatkinso/PGC/Afr-Eur-Ref/AFR_EUR.1kg_order.indivs.pops.txt -g /home/eatkinso/PGC/Afr-Eur-Ref/HapMapcomb_genmap_chr${i}.txt -n 5 -o PGC-SCZ.Bigdeli_merged.aa.rfmix.chr$i ;done



##and then run the Tractor script
#can launch this in pieces as the sections finish. Don't need fb results so actually can launch before those complete. But need to change to be pipes first. Also would be a lot faster seems like if split by chr first. Do this in tandem. 

for i in {1..1}; do python /home/eatkinso/PGC/Tractor/ExtractTracts.py --msp PGC-SCZ.Bigdeli_merged.aa.rfmix.chr$i --vcf PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr$i ;done
gsutil mv PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr1.*vcf gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schizophrenia/Tractor_output/




#combine these. Hapcount on Tractor VM
awk '{if ($1 ~ "#CHROM") print $0}' Bigdeli_VCF_header.txt > Bigdeli_VCF_header1.txt



cat Bigdeli_VCF_header1.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr1.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr2.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr3.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr4.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr5.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr6.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr7.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr8.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr9.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr10.anc0.dosage.txt  PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr11.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr12.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr13.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr14.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr15.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr16.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr17.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr18.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr19.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr20.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr21.anc0.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr22.anc0.dosage.txt | bgzip > PGC-SCZ.Bigdeli_merged.aa.anc0.dosage.txt.gz

cat Bigdeli_VCF_header1.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr1.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr2.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr3.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr4.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr5.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr6.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr7.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr8.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr9.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr10.anc1.dosage.txt  PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr11.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr12.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr13.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr14.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr15.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr16.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr17.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr18.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr19.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr20.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr21.anc1.dosage.txt PGC-SCZ.Bigdeli_merged2.aa.MAF1.chr22.anc1.dosage.txt | bgzip > PGC-SCZ.Bigdeli_merged.aa.anc1.dosage1.txt.gz



#then export again to the cloud
gsutil cp PGC-SCZ.Bigdeli_merged1.aa.anc0.dosage.txt.gz gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schizophrenia/Tractor_output/
gsutil cp PGC-SCZ.Bigdeli_merged1.aa.anc1.dosage.txt.gz gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schizophrenia/Tractor_output/

gsutil cp gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schizophrenia/Tractor_output/PGC-SCZ.Bigdeli_merged.aa.anc0.hapcount.txt.gz PGC-SCZ.Bigdeli_merged.aa.anc0.hapcount.txt.gz
gsutil cp PGC-SCZ.Bigdeli_merged.aa.anc1.hapcount.txt.gz gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schizophrenia/Tractor_output/

gsutil cp PGC-SCZ.Bigdeli_merged.aa.anc1.hapcount.txt.gz gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schizophrenia/Tractor_output/
gsutil cp PGC-SCZ.Bigdeli_merged.aa.anc0.hapcount1.txt.gz gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schizophrenia/Tractor_output/



#now can load into and run linreg in Hail
gsutil cp PGC-SCZ.Bigdeli_merged.aa.anc0.dosage.txt.gz gs://ukb-diverse-pops/AdmixedAfrEur/PGC-Schi
zophrenia/Tractor_output/

hailctl dataproc start ega --region us-central1  --packages gnomad

hailctl dataproc connect ega notebook

gcloud dataproc clusters update ega --num-preemptible-workers 10 --region us-central1

hailctl dataproc stop ega --region us-central1
