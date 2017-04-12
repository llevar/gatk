## Obtaining the datasets

Many can be obtained from the [GATK resource bundle](https://software.broadinstitute.org/gatk/download/bundle).

### Convert fasta to .2bit

```bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
gunzip human_g1k_v37.fasta.gz
faToTwoBit human_g1k_v37.fasta human_g1k_v37.2bit
```