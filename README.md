# nf-clair3-snp-calling

## Command

```bash
nextflow run main.nf -profile slurm --samplesheet /scratch/tweber/DATA/HGSVC_POOLS_PROJECT/DEMULTIPLEXING_POOLS/1000G_SNV_with_GT/NON_REFERENCED_SAMPLES/samplesheet_updated.csv --output_dir /scratch/tweber/DATA/HGSVC_POOLS_PROJECT/DEMULTIPLEXING_POOLS/1000G_SNV_with_GT/NON_REFERENCED_SAMPLES/clair3-results --model_name r941_prom_hac_g360+g422 --platform ont --reference /g/impont/ref/hg38.fa
```