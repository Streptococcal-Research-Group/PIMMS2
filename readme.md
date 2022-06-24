# PIMMS2

pimms2.py is a Python program to extract mutation/insertion sites/counts from sequencing data produced from a  [PIMMS](https://www.frontiersin.org/articles/10.3389/fmicb.2016.01645/full) experiment 


## Installation

Use the provided environment.yml file in conjunction with conda to set up a virtual environment for use with pimms2.py.

```bash
conda env create -f environment.yml
```
To start the environment:
```bash
conda activate venv_pimms2
```


## Usage

```python
python pimms2.py -h
usage: pimms2 [-h] [-v] {find_flank,bam_extract,table_merge,full_process} ...

description here

positional arguments:
  {find_flank,bam_extract,table_merge,full_process}
    find_flank          Mode: find read regions flanking the IS sequence by
                        mapping them to the target genome
    bam_extract         Mode: extract insertion site coordinates from sam file
    table_merge         Mode: merge two compatible PIMMS results tables (N.B:
                        this step does a simple table join and does not check
                        the data)
    full_process        Mode: find_flank + bam_extract

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

```

## Contributing
Open an issue or contact the developers to discuss requested changes.


## License
???