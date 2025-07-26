

~~~sh

dir=~/Bodymeta
cd $dir
##downlaod data
cat seq/fastq.txt|while read id
do
ascp -QT -l 300m -k1 -P33001  \
-i  ~/miniconda3/envs/download/etc/asperaweb_id_dsa.openssh   \
era-fasp@$id  seq/
done

##Import data
awk 'NR==1{print "SampleID\tforward-absolute-filepath\treverse-absolute-filepath"} \
  NR>1{print $1"\t$PWD/seq/"$1"_1.fastq.gz\t$PWD/seq/"$1"_2.fastq.gz"}' \
  result/metadata.txt > temp/manifest.txt
#cd ..
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $dir/temp/manifest.txt \
  --output-path $dir/temp/demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
qiime demux summarize  \
  --i-data $dir/temp/demux.qza \
  --o-visualization $dir/temp/demux.qzv
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs temp/demux.qza \
  --p-n-threads 10 \
  --p-max-ee-f 2 --p-max-ee-r 2 \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-trunc-len-f 250 --p-trunc-len-r 240 \
  --o-table temp/table.qza \
  --o-representative-sequences temp/rep-seqs.qza \
  --o-denoising-stats temp/denoising-stats.qza
qiime metadata tabulate \
  --m-input-file temp/denoising-stats.qza \
  --o-visualization temp/denoising-stats.qzv
qiime feature-table summarize \
  --i-table temp/table.qza \
  --o-visualization temp/table.qzv \
  --m-sample-metadata-file ./metadata.txt
qiime feature-table tabulate-seqs \
  --i-data temp/rep-seqs.qza \
  --o-visualization temp/rep-seqs.qzv
qiime feature-classifier classify-sklearn \
  --i-classifier data/classifier/classifier-V3V4.qza \
  --i-reads temp/rep-seqs.qza \
  --o-classification temp/taxonomy.qza
#Calculate the evolutionary tree
time qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences temp/rep-seqs.qza \
  --o-alignment temp/aligned-rep-seqs.qza \
  --o-masked-alignment temp/masked-aligned-rep-seqs.qza \
  --o-tree temp/unrooted-tree.qza \
  --o-rooted-tree temp/rooted-tree.qza
#Export the results of qiime2
qiime tools export --input-path temp/table.qza --output-path phyloseq/
qiime tools export --input-path temp/taxonomy.qza --output-path phyloseq/
qiime tools export --input-path temp/unrooted-tree.qza --output-path phyloseq/
qiime tools export --input-path temp/rep-seqs.qza --output-path phyloseq/
#The header of the taxonomy file must be converted to match proper format for merging with biom file
sed 's/Feature ID/#OTUID/' phyloseq/taxonomy.tsv | sed 's/Taxon/taxonomy/' | sed 's/Confidence/confidence/' > phyloseq/biom-taxonomy.tsv
rm phyloseq/taxonomy.tsv
biom add-metadata \
    -i phyloseq/feature-table.biom \
    -o phyloseq/taxa_table.biom \
    --observation-metadata-fp phyloseq/biom-taxonomy.tsv \
    --sc-separated taxonomy
#Functional prediction
time picrust2_pipeline.py -s phyloseq/dna-sequences.fasta \
    -i phyloseq/feature-table.biom \
    -o ./picrust2 \
    -p 10  
#rm -rf ../seq
#downstream analys
Rscript down_stream.R
~~~

