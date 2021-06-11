# enrichmenTE
Software for the analysis of LTR retrotransposon insertions in PCR-enriched NGS samples

# Usage
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