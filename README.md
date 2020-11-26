# Pragmatic Insertion Mutant Mapping System V2

## Getting Started

### Creating Conda Environment
    conda create --name pimms2 python=3.6
    conda activate pimms2
    conda install pip
  
### Install Dependancies  
    pip install -r requirements.txt
    
### Simple usage    

    conda activate pimms2  

    python ./pimms2_find_flanking.py find_flank -c pimms2.config --fasta <name>.fasta --mapper bwa --min 20 --max 50 --in_dir <name> --out_dir <name> --label test  

    python ./pimms2_sam_coords.py sam_extract -c pimms2.config --sam test/<name>_test_pimms2out_trim50_sub1_bwa.sam --label test --mismatch 0 --min_depth 1 --gff <name>.gff3  


### Options and Usage

    python ./pimms2_find_flanking.py find_flank -h  
    
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

    python ./pimms2_sam_coords.py sam_extract -h  
    
    Arguments:  
    -h, --help            show this help message and exit  
    -c --config           pimms2.config use parameters from config file  
    --sam pimms.sam/bam   sam/bam file of mapped IS flanking sequences  
    --label condition     text tag to add to results file  
    --mismatch float      fraction of permitted mismatches in mapped read ( 0 <= float < 0.2 [no filter]  
    --min_depth int       minimum read depth at insertion site >= int [2]  
    --gff genome.gff      GFF3 formatted file to use (note fasta sequence present in the file must be deleted before use)  
    --gff_extra 'x,y,z'   comma separated list of extra fields to include from the GFF3 annotation e.g. 'ID,translation,note'  
    --gff_force           override GFF/BAM seq id discrepancies e.g. use when the gff has a plasmid not present in the reference sequence or vice-versa  


