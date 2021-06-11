# enrichmenTE
Software for the analysis of LTR retrotransposon insertions in PCR-enriched NGS samples

# Table of Contents
* [Installation](#install)
* [Usage](#run)
* [Getting help](#help)
* [Citation](#citation)
* [License](#license)

# <a name="install"></a> Installation
## Use Conda to install software dependencies
enrichmenTE is written in python3 and is designed to run on a Linux operating system. Installation of software dependencies for enrichmenTE is automated by Conda, thus a working installation of Conda is required to install enrichmenTE. Conda can be installed via the Miniconda installer.

### Installing Miniconda (Python 3.X)
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc
conda init
```
- `conda init` requires you to close and open a new terminal before it take effect

### Update Conda
```
conda update conda
```

## Install software dependencies
After installing and updating Conda, you can now use conda to install dependencies and create running environment for enrichmenTE.
### Clone enrichmenTE Repository
```
git clone git@github.com:bergmanlab/enrichmenTE.git
cd enrichmenTE
```
### Create enrichmenTE Conda Environment
```
conda env create -f env/enrichmenTE.yml
```
- This installs all the software dependencies needed to run enrichmenTE into the enrichmenTE Conda environment.

### Activate enrichmenTE Conda Environment
```
conda activate enrichmenTE
```
- This adds the dependencies installed in the enrichmenTE conda environment to the environment PATH so that they can be used by the enrichmenTE scripts.
- This environment must always be activated prior to running any of the enrichmenTE scripts
- NOTE: Sometimes activating conda environments does not work via conda activate myenv when run through a script submitted to a queueing system, this can be fixed by activating the environment in the script as shown below
```
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate enrichmenTE
```
- For more on Conda: see the [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/index.html).


# <a name="run"></a> Usage
```
usage: enrichmenTE.py [-h] {detect,cluster} ...

Software to detect non-reference TEs from TE-NGS data enrichmenTE consists of two major steps: -
DETECT detects non-reference insertions for five TE families: 1731,297,copia,mdg1,roo - CLUSTER
cluster samples based on TE profiles of five TE families

positional arguments:
  {detect,cluster}  modes
    detect          Detect TE insertions using TE-NGS data
    cluster         Cluster samples based on TE profile

optional arguments:
  -h, --help        show this help message and exit
```

# <a name="help"></a> Getting help
Please use the [Github Issue page](https://github.com/bergmanlab/enrichmenTE/issues) if you have questions.

# <a name="citation"></a> Citation
To cite enrichmenTE in publications, please use https://www.biorxiv.org/content/10.1101/2021.04.24.441253v1

# <a name="license"></a> License
Copyright (c) 2020 Shunhua Han and Casey M. Bergman

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
