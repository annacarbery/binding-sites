conda activate esm_env
python src/predict_residues.py -t $1
conda activate pymol_env
python src/predict_centres.py -t $1
