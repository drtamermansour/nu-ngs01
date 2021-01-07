# Trimmomatic Tutorial

## Install the software

```
conda activate ngs1
conda install -c bioconda trimmomatic 
```

## Run Trimmomatic

### 1. Create practicing directory

```bash
mkdir -p ~/workdir/trimmed && cd ~/workdir/trimmed 
```

### 2. Set environmental variables

```bash
f1="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.pe.fq.gz"
f2="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.pe.fq.gz"

newf1="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.pe.trim.fq.gz"
newf2="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.pe.trim.fq.gz"

newf1U="$HOME/workdir/fqData/BD143_TGACCA_L005_R1_001.se.trim.fq.gz"
newf2U="$HOME/workdir/fqData/BD143_TGACCA_L005_R2_001.se.trim.fq.gz"

adap="$CONDA_PREFIX/share/trimmomatic-0.39-1/adapters"
```

### 3. Run trimmomatic

```bash
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:$adap/TruSeq3-PE-2.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36
```


Check the files in the adaptor folder (Source and more discussion at: https://www.biostars.org/p/323087/). You can see:
1. TruSeq2-PE.fa: Old libraries (usually data from the GAII), 
2. TruSeq3-PE.fa: libraries from the HiSeq or later machines until Nextera
3. TruSeq2-SE.fa and ruSeq3-SE.fa: as above but for SE libraries
4. TruSeq3-PE-2.fa: contains some additional sequences which find partial adapters in unusual location/orientation.
5. NexteraPE-PE.fa: as TruSeq3-PE-2.fa but for Nextera


Note: Specifying a trimlog file creates a log of all read trimmings, indicating the following details:
1. the read name
2. the surviving sequence length
3. the location of the first surviving base, aka. the amount trimmed from the start
4. the location of the last surviving base in the original read
5. the amount trimmed from the end

