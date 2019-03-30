
# Method One using our program GTF-Compare

### Create virtual evironment with conda
```
conda create -n ngs-gtf python=3.6 anaconda
source activate ngs-gtf
conda install -c conda-forge pypy3.5
wget https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py
```
### Install prerequisites

```
pypy3 -m pip install gffutils numpy tqdm 'intervaltree<3.0'
```

### Download required files
```
mkdir -p ~/workdir/gtf-compare/gtfs && cd ~/workdir/gtf-compare/gtfs
ln -s ~/workdir/hisat_align/ref_free.gtf .
ln -s ~/workdir/hisat_align/ref_sup.gtf .
mkdir -p ~/workdir/gtf-compare/method_one && cd ~/workdir/gtf-compare/method_one
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/comp.py
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/stat.py
```
### Run
```
source activate ngs-gtf
cd ~/workdir/gtf-compare/method_one
pypy3 comp.py -r ../gtfs/ref_sup.gtf ../gtfs/ref_free.gtf
pypy3 stat.py
```
###-----------------------------------------

# Method Two using GFFCompare

### Install gffcompare
```
source activate ngs-gtf
conda install gffcompare
```

### Run
```
mkdir -p ~/workdir/gtf-compare/method_two && cd ~/workdir/gtf-compare/method_two
gffcompare -r ../gtfs/ref_sup.gtf ../gtfs/ref_free.gtf
```
