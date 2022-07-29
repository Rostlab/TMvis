# AlphaFold 2 structures and TMbed predictions

**_TMvis_** requires two input files
1. AlphaFold 2 structures in `PDB` format (`tar`)
2. TMbed predictions (`txt`), see data format below

## Human proteome
To apply **_TMvis_** for the entire human proteome, run
```
chmod a+rx download_data.sh
./download_data.sh
```
from this directory. The script automatically creates a new folder, and downloads the data described below. This requires about 5GB storage capacity.

## Other use cases
To run **_TMvis_** for other proteomes, you can change `./download_data.sh` according to your needs.
- The available AlphaFold 2 structures can be found here [https://alphafold.ebi.ac.uk/download](https://alphafold.ebi.ac.uk/download). You must copy the download link into `./download_data.sh`
- TMbed predictions for all sequences in UniProtKB/Swiss-Prot (11-05-2022) are available [here](https://rostlab.org/public/tmbed/predictions/swissprot_110522_tmbed.tar.gz) (133MB)
- To run TMbed yourself, follow instructions [here](https://github.com/Rostlab/TMbed). This will help you create the text file (default TMbed settings) you need to run **_TMvis_**.

##
### TMbed data format 
- predicted residue classes are encoded by single letters
- 3-line format, every protein is represented by three lines: header, sequence, labels


`B`: Transmembrane beta strand (IN-->OUT orientation)\
`b`: Transmembrane beta strand (OUT-->IN orientation)\
`H`: Transmembrane alpha helix (IN-->OUT orientation)\
`h`: Transmembrane alpha helix (OUT-->IN orientation)\
`S`: Signal peptide\
`.`: Non-Transmembrane\

```
>7acg_A|P18895|ALGE_PSEAE
MNSSRSVNPRPSFAPRALSLAIALLLGAPAFAANSGEAPKNFGLDVKITGESENDRDLGTAPGGTLNDIGIDLRPWAFGQWGDWSAYFMGQAVAATDTIETDTLQSDTDDGNNSRNDGREPDKSYLAAREFWVDYAGLTAYPGEHLRFGRQRLREDSGQWQDTNIEALNWSFETTLLNAHAGVAQRFSEYRTDLDELAPEDKDRTHVFGDISTQWAPHHRIGVRIHHADDSGHLRRPGEEVDNLDKTYTGQLTWLGIEATGDAYNYRSSMPLNYWASATWLTGDRDNLTTTTVDDRRIATGKQSGDVNAFGVDLGLRWNIDEQWKAGVGYARGSGGGKDGEEQFQQTGLESNRSNFTGTRSRVHRFGEAFRGELSNLQAATLFGSWQLREDYDASLVYHKFWRVDDDSDIGTSGINAALQPGEKDIGQELDLVVTKYFKQGLLPASMSQYVDEPSALIRFRGGLFKPGDAYGPGTDSTMHRAFVDFIWRF
SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.........BBBBBBBBBB.................bbbbbbbbbbb.....BBBBBBBBBB...............................bbbbbbbbbb.........BBBBBB...............bbbbbbbb....BBBBBBBB....................bbbbbbbb......BBBBBBBB..........................bbbbbbbb..........BBBBBBBBBB............................bbbbbbbbbb.....BBBBBBBB..............................................bbbbbbbbb.....BBBBBBBB............................bbbbbbbbb..................BBBBBBBBBB...............bbbbbbbbb.
```
