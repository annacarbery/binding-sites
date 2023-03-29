import os
import json 
import matplotlib.pyplot as plt 

results = {str(num):{} for num in [1, 5, 10, 20, 30, 40]}
targets = json.load(open('binding-sites/data/test_set_A.json', 'r'))

pdb = []
pred = []

both_succeed = []
af2_fail = []
xtal_fail = []
both_fail = []


for num in results:
    af2 = {'0':[], '1':[], '2':[]}
    xtal = {'0':[], '1':[], '2':[]}
    for target in targets:
        if os.path.exists(f'binding-sites/predictions/{target}/success_af2_{num}.txt'):
            results_af2 = open(f'binding-sites/predictions/{target}/success_af2_{num}.txt', 'r').readlines()
            succ = False
            for rank in range(3):
                if len(results_af2) > rank and results_af2[rank].strip() == 'True':
                    succ = True
                af2[str(rank)].append(succ)

        if os.path.exists(f'binding-sites/predictions/{target}/success_xtal_{num}.txt'):
            results_xtal = open(f'binding-sites/predictions/{target}/success_xtal_{num}.txt', 'r').readlines()
            succ = False
            for rank in range(3):
                if len(results_xtal) > rank and results_xtal[rank].strip() == 'True':
                    succ = True
                xtal[str(rank)].append(succ)
        
        if num == '40':
            if af2['0'][-1] == True and xtal['0'][-1] == True:
                both_succeed.append(target)
            elif xtal['0'][-1] == True:
                af2_fail.append(target)
            elif af2['0'][-1] == True:
                xtal_fail.append(target)
            else:
                both_fail.append(target)



    results[str(num)]['af2'] = [af2[str(rank)].count(True)/len(af2[str(rank)]) for rank in range(3)]
    results[str(num)]['xtal'] = [xtal[str(rank)].count(True)/len(xtal[str(rank)]) for rank in range(3)]

    print(num)
    print(f'xtal success of {len(af2["0"])} targets:', results[str(num)]['xtal'])
    print(f'af2 success of {len(af2["0"])} targets:', results[str(num)]['af2'])
    
    pdb.append(results[str(num)]['xtal'][0])
    pred.append(results[str(num)]['af2'][0])

plt.plot([1, 5, 10, 20, 30, 40], pdb, label='PDB')
plt.plot([1, 5, 10, 20, 30, 40], pred, label='AF2')
plt.legend()
plt.ylabel('top rank success rate')
plt.xlabel('number of predictive models')
plt.tight_layout()
plt.savefig('binding-sites/plots/combination_success.png')

print(len(both_succeed), len(af2_fail), len(xtal_fail), len(both_fail))

lig_sim = json.load(open('HOLO4K/datafiles/dataset_similarity/ligand_similarity.json', 'r'))
pocket_sim = json.load(open('PocketShape/similarity_scores7.json', 'r'))

both_succeed_lig = []
af2_fail_lig = []
xtal_fail_lig = []
both_fail_lig = []

both_succeed_pock = []
af2_fail_pock = []
xtal_fail_pock = []
both_fail_pock = []

for targ in lig_sim:
    lig_sims = []
    lig_mcss = []

    for train in lig_sim[targ]:
        if len(lig_sim[targ][train]) > 0:
            lig_sims.append(max(lig_sim[targ][train]))
    
    best_sims = [pocket_sim[t] for t in pocket_sim if targ in t]+[0]

    if len(lig_sims) > 0:
        if targ in both_succeed:
            both_succeed_lig.append(max(lig_sims))
            both_succeed_pock.append(max(best_sims))
        elif targ in af2_fail:
            af2_fail_lig.append(max(lig_sims))
            af2_fail_pock.append(max(best_sims))
        elif targ in xtal_fail:
            xtal_fail_lig.append(max(lig_sims))
            xtal_fail_pock.append(max(best_sims))
        else:
            both_fail_lig.append(max(lig_sims))
            both_fail_pock.append(max(best_sims))



    


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
ax1.violinplot([both_succeed_lig, af2_fail_lig, xtal_fail_lig, both_fail_lig], showmeans=True)
ax1.set_xticks([1, 2, 3, 4])
ax1.set_xticklabels(['both succeed', 'af2 fail', 'pdb fail', 'both fail'])
ax1.set_ylabel('ligand similarity to training set')

ax2.violinplot([both_succeed_pock, af2_fail_pock, xtal_fail_pock, both_fail_pock], showmeans=True)
ax2.set_xticks([1, 2, 3, 4])
ax2.set_xticklabels(['both succeed', 'af2 fail', 'pdb fail', 'both fail'])
ax2.set_ylabel('pocket similarity to training set')

plt.tight_layout()
plt.savefig('binding-sites/plots/ligand_pocket_similarity.png')