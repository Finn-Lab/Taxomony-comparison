# Taxonomy-comparison
Repo for comparing reads taxonomy annotations between MG-RAST and MGnify
Comparison of taxonomies MG-RAST and MGnify

### RUNNING ==============================

1) MGnify <br>
Run: 
```
/hps/nobackup2/production/metagenomics/pipeline/tools-v4/miniconda2-4.0.5/bin/python  /hps/nobackup2/production/metagenomics/production-scripts/current/mgportal/analysis-pipeline/python/pipelineInitialisation.py -s /hps/nobackup2/production/metagenomics/production-scripts/current/mgportal/analysis-pipeline/python -f recalc -l 100 -y /hps/nobackup2/production/metagenomics/pipeline/tools/miniconda2-4.0.5/bin/python -o <outdir> -p <fasta> OR <1.fq,2.fq>
```
Necessary files: 
* taxonomy-results/SSU/SRR6367227_MERGED_FASTQ_SSU.fasta.mseq

2) MG-RAST <br>
For FASTA - run amplicon-fasta.workflow.cwl <br>
For FASTQ - run amplicon-fastq.workflow.cwl !!! Run with interleaved fastq file !!! <br>
Run: 
- Install from GitHub
- Create venv with mini-conda
- bsub -M 5000 -Is $SHELL
- activate vent
- add path of fasta/fastq file to yaml file
- run with singularity: 
```
cwltool --cachedir .cache --singularity --no-match-user amplicon-fastq.workflow.cwl amplicon-fastq.job.yaml
or 
cwltool --cachedir .cache --singularity --no-match-user amplicon-fasta.workflow.cwl amplicon-fasta.job.yaml 
```

Necessary files:  <br>
* amp_fq_test.440.cluster.rna97.mapping - list of clusters with all members
* amp_fq_test.440.cluster.rna97.fna - annotation of main members of cluster

### PARSING ==============================

Run python_comparison.py

### PLOT ==============================

Copy all lines after http://sankeymatic.com/build/ to web-visualiser


==============================
## Main ideas of comparison

The comparison is doing by each level and between reads that were annotated by both pipelines. The first step is to calculate the number of reads which were annotated further than super kingdom for MGnify and MG-RAST, calculate the number of reads that have their annotation stoped on sk level for MGnify and MG-RAST. Further look only on reads that were annotated deeper than sk. Repeat for kingdom, phylum and so on.

0) MG-RAST pipeline: <br>
Add to pipeline CWL file lines to output necessary files!
MG-RAST makes annotations for main members of cluster. All members in cluster have the same taxonomy as main member.

1) Kingdom <br>
Sometimes could be absent in taxonomic trees. It seems that MG-RAST skips this level and make annotation {super kingdom, phylum, class,…}. MGnify annotates “k__”. 
Solution: 
* add to all MG-RAST annotations prefixes “sk__, k__, and so on”
* add empty kingdom after super kingdom to MG-RAST annotation

2) MG-RAST (class) <br>
Sometimes MG-RAST duplicate class to further field, that do not have annotation. It is necessary to check names of fields excluding “(class)” from name. For example, NAME and NAME (CLASS) - this annotations are the same, but if we want to compare strings, they will be different.

3) MG-RAST Unclassified <br>
MG-RAST likes to annotate some fields as “unclassified came from …”. This annotations mean nothing - skip them. Be attentive with cases : sk__…;k__unclassified;p__unclassified;c_CLASSIFIED. This means that some fields absent in taxonomic tree. These cases must continue their comparison.

4) MG-RAST uncultured <br>
Skip all annotations with uncltured. 

5) Levenshtein distance <br>
Some field in MGnify and MG-RAST taxonomies could be different on one/two letters: c__Fusobacteriia", "c__Fusobacteria
It makes sense to calculate Levenshtein distance between lines. Lets say that annotations are different if LD > 3.

