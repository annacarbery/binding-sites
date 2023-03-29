from pymol import cmd
import json
import numpy as np
from sklearn.cluster import DBSCAN
import os
import pandas as pd
import time 
import argparse

def write_xyz(coords, target):
    """
    writes point cloud of pseudoatoms boxed around protein into xyz file
    """
    
    min_coords = [np.min(coords[:,0])-10, np.min(coords[:,1])-10, np.min(coords[:,2])-10]
    max_coords = [np.max(coords[:,0])+10, np.max(coords[:,1])+10, np.max(coords[:,2])+10]

    [x, y, z] = [np.arange(min_coords[i], max_coords[i], 1.5) for i in range(3)]

    with open(f'xyz/xyz_{target}.xyz', 'w') as w:
        w.write(f'{len(x)*len(y)*len(z)}\npoint\n')
        for a in x:
            for b in y:
                for c in z:
                    a_val = round(a, 3)
                    b_val = round(b, 3)
                    c_val = round(c, 3)
                    w.write(f'PS {a_val:.3f} {b_val:.3f} {c_val:.3f}\n')


def write_repeats(coords):
    with open('repeats.xyz', 'w') as w:
        w.write(f'{len(coords)}\npoint\n')
        for c in coords:
            w.write(f'PS {c[0]:.3f} {c[1]:.3f} {c[2]:.3f}\n')

def write_coords(coords, filename):
    with open(filename, 'w') as w:
        w.write(f'{len(coords)}\npoint\n')
        for c in coords:
            w.write(f'PS {c[0]:.3f} {c[1]:.3f} {c[2]:.3f}\n')


def get_final_cloud(pdb_xtal, pdb_af2, af2, target, final_preds, chain):
    """
    extracts points in point cloud that are close to protein residues that have been predicted to be druggable
    """
    cmd.load(pdb_xtal, 'complex')
    cmd.extract('hets', 'complex and HETATM')
    cmd.delete('hets')
    cmd.extract('chA', f'complex and chain {chain}')
    if af2 == True:
        cmd.delete('chA')
        cmd.load(pdb_af2, 'chA')

    coords = np.array(cmd.get_coords('chA')) 
    write_xyz(coords, target)
    assert os.path.exists(f'xyz/xyz_{target}.xyz')
    cmd.load(f'xyz/xyz_{target}.xyz','cloud1')

    cmd.extract('cloud2', 'cloud1 within 3 of prot') # don't want these
    cmd.extract('cloud3', 'cloud1 within 3 of chA')
    cmd.extract('cloud_bubble', 'cloud1 within 6 of chA') # don't want these

    for p in range(len(final_preds)):
        cmd.select(f'sele_{final_preds[p]}', f'chA and resi {final_preds[p]}') 
        if cmd.count_atoms(f'cloud_bubble within 4.5 of sele_{final_preds[p]}') > 0:
            cmd.create('cloud', f'cloud_bubble within 4.5 of sele_{final_preds[p]}')
            break

    for i in final_preds[p+1:]:

        try:
            cmd.select(f'sele_{i}', f'chA and resi {i}')
            cmd.select(f'cloud1_{i}', f'cloud_bubble within 4.5 of sele_{i}')
            if cmd.count_atoms(f'cloud1_{i}') > 0:
                cmd.copy_to('cloud', f'cloud1_{i}')
        except:
            raise


    coords = cmd.get_coords('cloud')
    coord_list = [','.join([str(n) for n in i]) for i in coords]
    repeats = list(sorted([i for i in coord_list if coord_list.count(i) > 2]))
    repeats = np.array([[float(x) for x in c.split(',')] for c in repeats])
    # write_repeats(repeats)

    return repeats


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target')
parser.add_argument('-num_models', '--models')
parser.add_argument('-af2', '--af2')
args = parser.parse_args()
target = args.target
num_models = int(args.models)
af2 = True if int(args.af2) == 1 else 0

pdb_xtal= f'HOLO4K/HOLO4k_data/orig_pdb/{target}.pdb'
strucs = [i for i in os.listdir('HOLO4K/datafiles/af2/aligned') if target in i]
pdb_af2 = f'HOLO4K/datafiles/af2/aligned/{strucs[0]}'
struc_type = 'xtal' if af2 == False else 'af2'

print(target)
[preds, chain] = json.load(open(f'binding-sites/predictions/{target}/predicted_residues_{struc_type}_{num_models}.json', 'r'))

site_success = []

try:
    cmd.reinitialize()
    final_coords = get_final_cloud(pdb_xtal, pdb_af2, af2, target, preds, chain)

    # cluster coordinates (that have been repeated by different residues) with a maximum distance threshold of 1.5A
    clustering = DBSCAN(eps=1.7, min_samples=2).fit(final_coords)

    # get size of each cluster for ranking
    site_counts = {}
    for site in set(clustering.labels_):
        if int(site) != -1:
            site_counts[site] = list(clustering.labels_).count(site)

    print('number of sites found:', len(site_counts))
    # determine centres of 3 largest potential binding sites and calculate whether ligand is near to binding site centre

    for site_rank in range(min([3,len(site_counts)])):
        biggest = max(site_counts, key=site_counts.get)

        # get coordinates of points in selected site, and calculate the centre by taking the mean of all points in the site
        main_site = np.array([final_coords[i] for i in range(len(final_coords)) if clustering.labels_[i] == biggest])
        write_coords(main_site, f'binding-sites/predictions/{target}/site_rank_{site_rank}_{struc_type}_{num_models}.xyz')
        centre = np.array([np.mean(main_site[:,0]), np.mean(main_site[:,1]), np.mean(main_site[:,2])])

        with open(f'binding-sites/predictions/{target}/centre_rank_{site_rank}_{struc_type}_{num_models}.xyz', 'w') as w:
            w.write(f'1\npoint\n')
            w.write(f'PS {centre[0]:.3f} {centre[1]:.3f} {centre[2]:.3f}\n')


        del site_counts[biggest]

        cmd.reinitialize()
        cmd.load(f'binding-sites/predictions/{target}/centre_rank_{site_rank}_{struc_type}_{num_models}.xyz', 'centre')
        cmd.load(pdb_xtal, 'complex')
        lig_names = [i.split('_')[0] for i in open(f'/data/xchem-fragalysis/tyt15771/projects/frag/HOLO4K/HOLO4k_data/data/{target}.ligand', 'r').readlines() if len(i.strip()) > 2]
        print('ligands:', list(set(lig_names)))
        cmd.extract('lig', f'complex and (resn {" or resn ".join(lig_names)})')
        success =  cmd.count_atoms('lig within 4 of centre') > 0
        print(site_rank, biggest, centre, success)
        site_success.append(str(success))
except:
    print(target, num_models, struc_type, 'error')

with open(f'binding-sites/predictions/{target}/success_{struc_type}_{num_models}.txt', 'w') as w:
    w.write('\n'.join(site_success))