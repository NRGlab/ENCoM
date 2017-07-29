from modeller import *
from modeller.automodel import *

import sys

prot = sys.argv[1]
prot_file = sys.argv[2]
env = environ()
env.io.atom_files_directory = ['.', '../atom_files']
aln = alignment(env)
env.io.hetatm = True
mdl = model(env, file=prot)
#aln.append_model(mdl, align_codes=prot, atom_files=prot_file)
#aln.append_model(mdl, align_codes='tobuild', atom_files=prot_file)
#aln.align()
#aln.write(file='all_aligne.pir', alignment_format='PIR')
a = automodel(env, alnfile='all_aligne.pir',
    knowns=prot, sequence='tobuild',
     assess_methods=assess.GA341
    )
    
    #assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = 1
a.final_malign3d = True
#a.very_fast()
a.make()

