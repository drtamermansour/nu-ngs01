# FASTQC Tutorial


## Download Fastq files
```
mkdir -p ~/workdir/fqData && cd ~/workdir/fqData
wget https://de.cyverse.org/dl/d/3CE425D7-ECDE-46B8-AB7F-FAF07048AD42/samples.tar.gz
tar xvzf samples.tar.gz
```

## Install the software
```
source activate ngs1
conda install -c bioconda fastqc 
conda install -c bioconda multiqc 
```

## Explore the help of FASTQC
```
mkdir -p ~/workdir/FASTQC_tut && cd ~/workdir/FASTQC_tut
fastqc -h > fastqc_hlp
less fastqc_hlp
```

## Run the FASTQC for each read end
```
for f in ~/workdir/fqData/*.fq.gz;do fastqc -t 1 -f fastq -noextract $f;done
```

## Merge the output reports into one super report
```
mv ../fqData/*html ./
mv ../fqData/*zip ./
multiqc -z -o . .
```

## View the results

1. Open jupyter
2. Browse to ~/workdir/FASTQC_tut
3. Double click on the file **multiqc_report.html**
4. On the top left of the multiqc_report.html tab header, click on **Trust HTML**
