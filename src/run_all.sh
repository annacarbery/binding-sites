source /data/xchem-fragalysis/tyt15771/miniconda3/bin/activate esm
cd /data/xchem-fragalysis/tyt15771/projects/frag/
python binding-sites/src/predict_residues.py -t $1 -num_models $2 -af2 0
python binding-sites/src/predict_residues.py -t $1 -num_models $2 -af2 1
source /data/xchem-fragalysis/tyt15771/miniconda3/bin/activate pymol
python binding-sites/src/predict_centres.py -t $1 -num_models $2 -af2 0
python binding-sites/src/predict_centres.py -t $1 -num_models $2 -af2 1
