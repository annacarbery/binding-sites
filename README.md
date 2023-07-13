# IF-SitePred

IF-SitePred is a method for predicting ligand-binding sites on protein structures. It first generates an embedding for each residue of the protein using the [ESM-IF1 (inverse folding)](https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding) model, then performs point cloud clustering to identify binding site centres. 

## Installation

Follow the following steps to prepare a virtual environment for IF-SitePred:

1. Install esm
```
pip install fair-esm
```
2. Install [torch](https://pytorch.org/) (system-dependent)
3. Install remaining dependencies
```
pip install scipy
pip install torch-geometric
pip install torch_scatter
pip install biotite
pip install lightgbm
```
4. Clone repository
```
git clone https://github.com/annacarbery/binding-sites
```

## Command line use

## Citation

In progress...
