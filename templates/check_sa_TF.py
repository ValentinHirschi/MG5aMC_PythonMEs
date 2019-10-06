#!/usr/bin/env python2

# This template is specialised for the TF output

import random
import os

from processes.all_processes import *
from model.parameters import ModelParameters
from phase_space_generator.flat_phase_space_generator import FlatInvertiblePhasespace 

class Colour:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

module_name = os.path.basename(os.path.dirname(os.path.realpath( __file__ )))

all_process_classes = [{all_process_classes}]

# For now, the feature of specifying an SLHA param card to initialise
# the value of independent parameters is not supported yet.
active_model = ModelParameters(None)

# Center of mass of the collision in GeV
E_cm = 14000.

print("")
print("The module '%s' contains %d processes"%(module_name, len(all_process_classes)))
print("")
print(str(active_model))
print("")

for process_class in all_process_classes:
    
    print(">>> Running process %s%s%s"%(Colour.BLUE,process_class.__name__,Colour.END))

    # Generate a random PS point for this process
    process = process_class()
    external_masses = process.get_external_masses(active_model)

    # Ensure that E_cm offers enough twice as much energy as necessary 
    # to produce the final states
    this_process_E_cm = max( E_cm, sum(external_masses[1])*2. )

    ps_generator = FlatInvertiblePhasespace(
        external_masses[0], external_masses[1],
        beam_Es = (this_process_E_cm/2.,this_process_E_cm/2.),
        # We do not consider PDF for this standalone check
        beam_types=(0,0)
    )

    # Generate some random variables
    random_variables = [random.random() for _ in range(ps_generator.nDimPhaseSpace())]

    PS_point, jacobian = ps_generator.generateKinematics(this_process_E_cm, random_variables)
    
    print("> PS point:")
    print(PS_point)
    print("> Matrix element evaluation : %s%.16e%s"%(Colour.GREEN,process.smatrix(PS_point, active_model),Colour.END))
    print("")

