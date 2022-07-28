# Pragmatic Insertion Mutant Mapping System V2

PIMMS2 is a Python program to extract mutation/insertion sites/counts from sequencing data produced from a  [PIMMS](https://www.frontiersin.org/articles/10.3389/fmicb.2016.01645/full) experiment 

## Getting Started

### Creating Conda Environment Linux and MacOS x86

Use the provided environment.yml file in conjunction with conda to set up a virtual environment for use with pimms2.py.

```bash
conda env create -f environment.yml
```
To start the environment:
```bash
conda activate venv_pimms2

### Creating a Conda environment MacOS M1

CONDA_SUBDIR=osx-64 conda create -n pimms2 python=3.6

conda actiavte pimms2

conda config --env --set subdir osx-64

pip install -r requirements.txt

    
### Simple usage    

	All parameters can be set in the pimms2.config file or you can specify them on the command line

    $ conda activate pimms2  

    $ python ./pimms2.py <module> -h

    $ python ./pimms2.py find_flank -c pimms2.config --fasta <name>.fasta --mapper bwa --min 20 --max 50 --in_dir <name> --out_dir <name> --label test  

    $ python ./pimms2.py bam_extract -c pimms2.config --bam test/<name>_test_pimms2out_trim50_sub1_bwa.bam --label test --mismatch 0 --min_depth 1 --gff <name>.gff3 
    
    $ python ./pimms2.py full_process -c pimms2.config 


### Options and Usage

    python ./pimms2.py find_flank -h  
    
    Arguments:  
    -h, --help            show this help message and exit  
    -c --config           pimms2.config use parameters from config file  
    --nano                override with settings more suitable for nanopore  
    --cas9                override with motif matching method settings for nanopore cas9  
    --fasta               fasta file for reference genome  
    --nomap               do not run mapping step  
    --mapper              minimap2 or bwa   
    --rmfiles             remove intermediate files   
    --lev                 use Levenshtein distance of 1  
    --sub SUB             number of permitted base substitutions in motif match [1]  
    --insert INSERT       number of permitted base insertions in motif match [0]  
    --del DELETION        number of permitted base insertions in motif match [0]  
    --in_dir IN_DIR       directory containing input fastq files (assumed to match '*fq.gz' or '*.fastq')  
    --fwdrev FWDREV       text substring to uniquely identify illumina fwd/rev paired fastq files ['_R1_,_R2_']  
    --out_dir DIR         directory to contain fastq files results  
    --cpus CPUS           number of processors to use [(os.cpu_count() / 2)]  
    --max MAX             clip results to this length [illumina:60/nano:100]  
    --min MIN             minimum read length [illumina:60/nano:100]  
    --motif1 MOTIF1       IS end reference motif1 [e.g. pGh9:TCAGAAAACTTTGCAACAGAACC]  
    --motif2 MOTIF2       IS end reference motif2 [e.g. pGh9:GGTTCTGTTGCAAAGTTTAAAAA]  
    --label condition     Text tag to add to results file  

    python ./pimms2.py bam_extract -h  
    
    Arguments:  
    -h, --help            show this help message and exit  
    -c --config           pimms2.config use parameters from config file  
    --bam pimms.bam/sam   bam/sam file of mapped IS flanking sequences  
    --label condition     text tag to add to results file  
    --mismatch float      fraction of permitted mismatches in mapped read ( 0 <= float < 0.2 [no filter]  
    --min_depth int       minimum read depth at insertion site >= int [2]  
    --gff genome.gff      GFF3 formatted file to use (note fasta sequence present in the file must be deleted before use)  
    --gff_extra 'x,y,z'   comma separated list of extra fields to include from the GFF3 annotation e.g. 'ID,translation,note'  
    --gff_force           override GFF/BAM seq id discrepancies e.g. use when the gff has a plasmid not present in the reference sequence or vice-versa  
	--out_fmt			  xlsx, tsv, csv

