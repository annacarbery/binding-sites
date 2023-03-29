import os
import json
import esm
import torch
import pickle
import numpy as np 
import argparse 


def get_chains(pdb_path):
    chains = [i[21:22] for i in open(pdb_path, 'r').readlines() if i.startswith('ATOM')]
    return list(set(chains))

def get_index_coords(pdb_path, chain):
    coords_resi = {}
    atoms = [i for i in open(pdb_path, 'r').readlines() if i.startswith('ATOM') and i[21:22] == chain]
    for atom in atoms:
        coords = [atom[30:38].strip(), atom[38:46].strip(), atom[46:55].strip()]
        coords_resi[','.join(coords)] = int(atom[22:26].strip())
    return coords_resi


def load_models(num_models):
    models = []
    for i in range(num_models):
        models.append(pickle.load(open(f'HOLO4K/models/lgbm_2803/lgbm_if1_{i}.pkl', 'rb')))
    return models


model, alphabet = esm.pretrained.load_model_and_alphabet_local('/data/xchem-fragalysis/tyt15771/projects/esm/esm_if1_gvp4_t16_142M_UR50.pt')
model.eval()

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target')
parser.add_argument('-num_models', '--models')
parser.add_argument('-af2', '--af2')
args = parser.parse_args()
target = args.target
num_models = int(args.models)
af2 = True if int(args.af2) == 1 else 0

if not os.path.isdir(f'binding-sites/predictions/{target}'):
    os.mkdir(f'binding-sites/predictions/{target}')

models = load_models(num_models)
preds = {}

if af2 == False:
    pdb = f'HOLO4K/HOLO4k_data/orig_pdb/{target}.pdb'
else:
    strucs = [i for i in os.listdir('HOLO4K/datafiles/af2/aligned') if target in i]
    pdb = f'HOLO4K/datafiles/af2/aligned/{strucs[0]}'

chains = sorted(get_chains(pdb))
preds[pdb] = {}

for chain in chains[:1]:
    print(chain)	

    structure = esm.inverse_folding.util.load_structure(pdb, chain)
    coords_resi = get_index_coords(pdb, chain)
    coords, seq = esm.inverse_folding.util.extract_coords_from_structure(structure)	
    if1 = esm.inverse_folding.util.get_encoder_output(model, alphabet, coords).detach().numpy()
    
    pred = []
    for model_flaml in models:
        p = np.rint(model_flaml.predict(if1))
        pred.append(p)

    pred_intersection = []
    model_num = len(models)
    for i in range(len(pred[0])):
        pred_intersection.append(1 if [pred[a][i] for a in range(model_num)].count(1) == model_num else 0)
    

    res_nums = []
    for i in range(len(coords)):
        coord_str = ','.join([str("{:.3f}".format(x)) for x in coords[i][0]])
        if pred_intersection[i] == 1:
            try:
                res_nums.append(str(coords_resi[coord_str]))
            except:
                pass
    preds[pdb][chain] = res_nums
    print(pdb, chain, res_nums, len(res_nums))
    if af2 == False:
        json.dump([res_nums, chain], open(f'binding-sites/predictions/{target}/predicted_residues_xtal_{num_models}.json', 'w'))
    else:
        json.dump([res_nums, chain], open(f'binding-sites/predictions/{target}/predicted_residues_af2_{num_models}.json', 'w'))