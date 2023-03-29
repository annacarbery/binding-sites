import json
import os 

test_targs = json.load(open('binding-sites/data/test_set_A.json', 'r'))
print('Test targets:', len(test_targs))

for num in [5, 10, 30]:
    for targ in test_targs:
        bash_path = f'/data/xchem-fragalysis/tyt15771/projects/frag/binding-sites/submissions/{targ}_{num}.sh'
        bash_template = open('binding-sites/src/run_all.sh', 'r').read()
        new_bash = bash_template.replace('$1', targ)
        new_bash = new_bash.replace('$2', str(num))
        with open(bash_path, 'w') as w:
            w.write(new_bash)

        submission_path = f'binding-sites/submissions/{targ}_{num}'
        submission_template = open('binding-sites/submissions/template', 'r').read()
        new_submission = submission_template.replace('BASH_PATH', bash_path)
        new_submission = new_submission.replace('TARGET', targ)
        new_submission = new_submission.replace('MODEL', str(num))
        with open(submission_path, 'w') as w:
            w.write(new_submission)

        os.system(f'condor_submit {submission_path}')
