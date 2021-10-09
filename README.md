# genome-comparator

Utility for [Genome browser](https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr8%3A128745680%2D128755674&hgsid=1181507859_VKASGgx2pB8RTrhe9OpaCVtWAsOS)
data comparison

## Environment preparation
`python3 -m pip install -r requirements.txt`

## Data preparation
In order to compare data, you need to:
  1. Download QmRLFS finder analysis data via Genome Browser API
  - e.g. with curl `curl -L 'https://api.genome.ucsc.edu/getData/track?genome=hg19;chrom=chr8;start=128745680;end=128755674;track=hub_124277_RLFS' > RLFS.txt`
  2. Download input DNA sequence for R-loop tracker analysis
  - e.g. with curl `curl -L 'https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr8;start=128745680;end=128755674;' > Neat1_chr8:128745680-128755674`
  3. Clean input sequence data (original response from Genome browser contains metadata aswell)
  - using `util.py` module: `python3 util.py` (you have to change the path of a input file)
  4. Run R-loop tracker analysis on [bioinformatics server](https://bioinformatics.ibp.cz) (follow [help instructions](https://bioinformatics.ibp.cz/#/help/rloopr))
  5. Download `bedgrahph` format, **change paths in the script** according to your data and run `genome_comparator.py`

This process will be automatized in the next version, you will only have to provide following params for the analysis: `chrom=chr8;start=128745680;end=128755674`
