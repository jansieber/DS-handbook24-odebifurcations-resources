# DS-handbook24-odebifurcations-resources
This repository contains computational scripts to reproduce the data for figures in

Harry Dankowicz \& Jan Sieber: Computational Bifurcation Analysis. Chapter in: Handbook on Nonlinear Dynamics. Volume 2 Numerical Methods. Edited by Vincent Acary and Oleg Gendelman.

## Folders:
* `coco_r3316`: snapshot of coco used for the runs. For newest stable version: (https://sourceforge.net/projects/cocotools/)
* `cstr:` computations for plots in Sections 2 and 3 on continuous stirred tank reactor
* `brus`: computations for section 4 on network of chemical oscillators
* `figures`: destinatino folder for `exportgraphics`

## Instructions

Computations require Matlab

1. Enter folder `coco_r3316`.
2. Execute `startup` script to set paths for `coco`.
3. Enter `cstr`.
4. Execute script `cstr_computations` to perform all computations. Results will be stored in subfolder `cstr/data`.
5. Execute `cstr_plot1`, `cstr_plot2` to regenerate figures from section 2 and 3.
6. Enter `brus`.
7. Execute scripts `branchoff_hopf`, `run_sympo`, `run_hopf` to perform all computations. Results will be stored in subfolder `brus/data`.
8. Execute scripts `branchoff_hopf_plot`, `sympo_plot`, `sbhopf_plot` to regenerate figures from section 4. 
