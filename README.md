# cfMEX_tools
cell-free DNA multi-feature extraction tools
  - [x] Bam to Fragment file
  - [x] Fragment size
  - [x] End motif
  - [x] Short per Long ratio
  - [ ] Nucleosome footprint
        
*uncheck boxes are currently under development*

# Installation
Check python version >= 3.12.0
```
python -V
```
Important modules
```
pip install numpy pandas pysam biopython
```

# Setting
You can change Reference Genome path align with your local
by changing in `cfMEX_tools/scr/__init__.py`

```
#######################################################################################
# SETTING

Reference_Genome = '/mnt/sas/ref/hg38/v0/Homo_sapiens_assembly38.fasta'

#######################################################################################
```
