# AI Usage Documentation


## Summary of Usage
Claude Sonnet 4.5 was used throughout this assignment for code generation and result interpretation. 

---

## Detailed Usage by Section

### 1. Reference Genome Download
  - Commands for downloading chr10 from UCSC/Ensembl
  - wget/curl syntax
  - FASTA file handling

### 2. Alignment (minimap2)
  - minimap2 parameter selection for Illumina vs PacBio data
  - Proper preset flags (-ax sr vs -ax map-pb)
  - SAM to BAM conversion and indexing commands
  

### 3. Variant Calling (bcftools)
  - bcftools mpileup and call command syntax
  - Parameter selection 
  - VCF filtering and normalization
  

### 4. Phasing (HapCUT2/WhatsHap)
  - Installation instructions for phasing tools
  - Command syntax for HapCUT2 and WhatsHap
  - Converting HapCUT2 block format to phased VCF
  - Troubleshooting phasing errors

### 5. Variant Analysis & Comparison
  - bcftools isec commands for VCF comparison
  - Extracting shared vs unique variants
  - IGV automation scripting
  - Interpreting IGV screenshots


### 6. Star-Allele Determination
  - Understanding PharmVar database structure
 
