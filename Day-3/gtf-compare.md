### Create virtual evironment with conda
```
conda create -n ngs-gtf python=3.6 anaconda
```
### Install prerequisites
```
source activate ngs-gtf
pip install gffutils
pip install numpy
pip install tqdm
pip install 'intervaltree<3.0'
```
### Download required files
```
mkdir ~/gtf-compare && cd ~/gtf-compare
mkdir databases && cd databases
wget https://transfer.sh/HjxD/databases.tar.xz
wget https://transfer.sh/QeKeX/gtfs.tar.xz
tar xvf databases.tar.xz
cd ..
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/comp.py
wget https://raw.githubusercontent.com/abdelrahmanMA/gtf-compare/master/code/stat.py
```
### Run
```
source activate ngs-gtf
cd ~/gtf-compare
python comp.py -r ./databases/gencode.v27.primary_assembly.annotation.gtf.db ./databases/ribo_tumor_scallop_merged.gtf.db
python comp.py -r ./databases/gencode.v27.primary_assembly.annotation.gtf.db ./databases/poly_tumor_scallop_merged.gtf.db
python stat.py
```
