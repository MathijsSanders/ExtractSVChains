# ExtractSVChains
Extract fold-back inversion and fold-back inversion with small templated insertion chains from BRASS output

## How do I run it?

First, BRASS has been run on your sample of interest and matched control for determining somatic structural variants (SV). Second, it is preferred that the SVs are filtered further by using AnnotateBRASS. Finally, the resultant BED file is used as input for this algorithm to determine the number of classical inversion events, fold-back inversions and fold-back-inversions supplemented with small templated insertions.

### The recommended way

The pre-compiled JAR file is included with the repository, but in case the package needs to be recompiled, please run:

```bash
mvn package clean
```

The following command adds additional statistics to the BRASS output for filtering purposes:

```bash
java -Xmx15G -jar extractChains.jar --input-bam-file input_bam_file  --input-brass-bam-file BRASS_bam_file --brass-bed-file filtered_bed_file --output-file output_file --reference reference_genome_FA_file --decoys comma_separated_names_of_decoys > additional_info_file
```

- --input-bam-file*: Input BAM file.
- --input-brass-bam-file*: Input BAM file produced by BRASS.
- --brass-bed-file*: Input BRASS BED output file.
- --output-file*: Output file for writing statistics.
- --reference*: Reference genome in FASTA format.
- --decoys: Comma-separated list of decoy names (e.g., hs37d5).    
- --width-mini-extract: Minimum window for extracting reads (mutation_position +- width, default: 20).
- --width-maxi-extract: Maximum window for extracting reads (mutation_position +- width, default: 1000).
- --width-minor-extract: Minor larger sized window for extracting reads (mutation_position +- width, default: 2000).
- --width-major-extract: Major window for extracting reads (mutation_position +- width, default: 20000).
- --help, -help: Get usage information
- --version, -version: Get current version
- \* Required.

*Dependencies*
- Maven version 3+ (For compiling only).
- Java JDK 1.8+
