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
pypy3 -m pip install gffutils
pypy3 -m pip install numpy
pypy3 -m pip install tqdm
pypy3 -m pip install 'intervaltree<3.0'
```
### Download required files
```
mkdir ~/gtf-compare && cd ~/gtf-compare
mkdir databases && cd databases
wget https://transfer.sh/HjxD/databases.tar.xz
wget https://transfer.sh/QeKeX/gtfs.tar.xz
tar -xvf databases.tar.xz
cd ..
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/comp.py
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/stat.py
```
### Run
```
source activate ngs-gtf
cd ~/gtf-compare
pypy3 comp.py -r ./databases/gencode.v27.primary_assembly.annotation.gtf.db ./databases/ribo_tumor_scallop_merged.gtf.db
pypy3 comp.py -r ./databases/gencode.v27.primary_assembly.annotation.gtf.db ./databases/poly_tumor_scallop_merged.gtf.db
pypy3 stat.py
```
