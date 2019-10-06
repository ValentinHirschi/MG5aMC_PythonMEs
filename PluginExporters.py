#####################################################
#                                                   #
#  Source file of the MG5aMC_PythonMEs plugin.      #
#                                                   #
#         author: Valentin Hirschi                  #
#                                                   #
#####################################################

import os
import stat
import logging
import itertools
from pprint import pformat
from math import fmod
import aloha
import shutil

import aloha.create_aloha as create_aloha

plugin_path = os.path.dirname(os.path.realpath( __file__ ))

from madgraph import MadGraph5Error, InvalidCmd, MG5DIR

import madgraph.iolibs.export_python as export_python
import madgraph.iolibs.export_cpp as export_cpp
import madgraph.core.base_objects as base_objects
from madgraph.iolibs.file_writers import PythonWriter

import madgraph.iolibs.helas_call_writers as helas_call_writers

import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.file_writers as writers
import madgraph.various.misc as misc
from madgraph.iolibs.files import cp, ln, mv

import madgraph.iolibs.helas_call_writers as helas_call_writers

logger = logging.getLogger('MG5aMC_PythonMEs.MEExporter')

pjoin = os.path.join 


class UFOModelConverterPython(export_cpp.UFOModelConverterCPP):

    # Static variables (for inheritance)
    output_name = 'Python Standalone'
    namespace = 'Python'

    # Dictionary from Python type to C++ type
    type_dict = {"real": "float",
                 "complex": "complex"}

    # Regular expressions for cleaning of lines from Aloha files
    compiler_option_re = None 
    namespace_re = None

    slha_to_depend = {('SMINPUTS', (3,)): ('aS',),
                      ('SMINPUTS', (1,)): ('aEM',)}

    # Template files to use
    include_dir = '.'
    cc_file_dir = '.'
    param_template_h = None
    param_template_cc = None
    aloha_template_h = None
    aloha_template_cc = None

    param_template_py = 'model_template.py'

    copy_include_files = []
    copy_cc_files = []

    def __init__(self, model, output_path, wanted_lorentz = [],
                 wanted_couplings = [], replace_dict={}):
        """ initialization of the objects """

        self.model = model
        self.model_name = export_cpp.ProcessExporterCPP.get_model_name(model['name'])
        self.aloha_model = create_aloha.AbstractALOHAModel(self.model_name)


        self.dir_path = output_path
        self.default_replace_dict = dict(replace_dict)
        # List of needed ALOHA routines
        self.wanted_lorentz = wanted_lorentz
        self.wanted_couplings = wanted_couplings

        # For dependent couplings, only want to update the ones
        # actually used in each process. For other couplings and
        # parameters, just need a list of all.
        self.coups_dep = {}    # name -> base_objects.ModelVariable
        self.coups_indep = []  # base_objects.ModelVariable
        self.params_dep = []   # base_objects.ModelVariable
        self.params_indep = [] # base_objects.ModelVariable
        self.p_to_cpp = None 

        # Prepare parameters and couplings for writeout in C++
        self.prepare_parameters()
        self.prepare_couplings(wanted_couplings)

    def write_files(self):
        """Create all necessary files"""

        os.makedirs(self.dir_path)
        open(pjoin(self.dir_path,'__init__.py'),'w')

        # Write Helas Routines
        self.write_aloha_routines()

        # Write parameter (and coupling) class files
        self.write_parameter_class_files()

    # Routines for preparing parameters and couplings from the model

    def prepare_parameters(self):
        """Extract the parameters from the model, and store them in
        the two lists params_indep and params_dep"""

        # Keep only dependences on alphaS, to save time in execution
        keys = self.model['parameters'].keys()
        keys.sort(key=len)
        params_ext = []
        for key in keys:
            if key == ('external',):
                params_ext += [p for p in self.model['parameters'][key] if p.name]
            elif 'aS' in key:
                for p in self.model['parameters'][key]:
                    self.params_dep.append(base_objects.ModelVariable(p.name,
                                              p.name + " = " + p.expr,
                                              p.type,
                                              p.depend))
            else:
                for p in self.model['parameters'][key]:
                    if p.name == 'ZERO':
                        continue
                    self.params_indep.append(base_objects.ModelVariable(p.name,
                                              p.name + " = " + p.expr,
                                              p.type,
                                              p.depend))

        # For external parameters, want to read off the SLHA block code
        while params_ext:
            param = params_ext.pop(0)
            # Read value from the slha variable
            expression = ""
            assert param.value.imag == 0
            if len(param.lhacode) == 1:
                expression = "%s = slha.get_block_entry(\"%s\", %d, %e);" % \
                             (param.name, param.lhablock.lower(),
                              param.lhacode[0], param.value.real)
            elif len(param.lhacode) == 2:
                expression = "indices[0] = %d;\nindices[1] = %d;\n" % \
                             (param.lhacode[0], param.lhacode[1])
                expression += "%s = slha.get_block_entry(\"%s\", indices, %e);" \
                              % (param.name, param.lhablock.lower(), param.value.real)
            else:
                raise MadGraph5Error("Only support for SLHA blocks with 1 or 2 indices")
            self.params_indep.insert(0,
                                   base_objects.ModelVariable(param.name,
                                                   expression,
                                                              'real'))
            
    def prepare_couplings(self, wanted_couplings = []):
        """Extract the couplings from the model, and store them in
        the two lists coups_indep and coups_dep"""

        # Keep only dependences on alphaS, to save time in execution
        keys = self.model['couplings'].keys()
        keys.sort(key=len)
        for key, coup_list in self.model['couplings'].items():
            if "aS" in key:
                for c in coup_list:
                    if not wanted_couplings or c.name in wanted_couplings:
                        self.coups_dep[c.name] = base_objects.ModelVariable(\
                                                                   c.name,
                                                                   c.expr,
                                                                   c.type,
                                                                   c.depend)
            else:
                for c in coup_list:
                    if not wanted_couplings or c.name in wanted_couplings:
                        self.coups_indep.append(base_objects.ModelVariable(\
                                                                   c.name,
                                                                   c.expr,
                                                                   c.type,
                                                                   c.depend))

        # Convert coupling expressions from Python to C++
        for coup in self.coups_dep.values() + self.coups_indep:
            coup.expr = coup.name + " = " + coup.expr

    # Routines for writing the parameter files

    def write_parameter_class_files(self):
        """Generate the python model parameters 
        which have the parameters and couplings for the model."""

        parameters_file_path = pjoin(self.dir_path,'parameters.py')
        parameters_file = PythonWriter(parameters_file_path)


        file_py = self.generate_parameters_class_files()

        # Write the files
        parameters_file.write(file_py)
        parameters_file.close()

        logger.info("Created parameters python file %s."%parameters_file_path)


    def generate_parameters_class_files(self, indent=8):
        """Create the content of the Parameters_model.h and .cc files"""

        replace_dict = self.default_replace_dict

        replace_dict['info_lines'] = export_cpp.get_mg5_info_lines()
        replace_dict['model_name'] = self.model_name

        replace_dict['independent_parameters'] = \
                                   "%s# Model parameters independent of aS\n"%(' '*indent) + \
                                   self.write_parameters(self.params_indep)
        replace_dict['independent_couplings'] = \
                                   "%s# Model parameters dependent on aS\n"%(' '*indent) + \
                                   self.write_parameters(self.params_dep)
        replace_dict['dependent_parameters'] = \
                                   "%s# Model couplings independent of aS\n"%(' '*indent) + \
                                   self.write_parameters(self.coups_indep)
        replace_dict['dependent_couplings'] = \
                                   "%s# Model couplings dependent on aS\n"%(' '*indent) + \
                                   self.write_parameters(self.coups_dep.values())

        replace_dict['set_independent_parameters'] = \
                               self.write_set_parameters(self.params_indep)
        replace_dict['set_independent_couplings'] = \
                               self.write_set_parameters(self.coups_indep)
        replace_dict['set_dependent_parameters'] = \
                               self.write_set_parameters(self.params_dep)
        replace_dict['set_dependent_couplings'] = \
                               self.write_set_parameters(self.coups_dep.values())

        replace_dict['print_independent_parameters'] = \
                               self.write_print_parameters(self.params_indep)
        replace_dict['print_independent_couplings'] = \
                               self.write_print_parameters(self.coups_indep)
        replace_dict['print_dependent_parameters'] = \
                               self.write_print_parameters(self.params_dep)
        replace_dict['print_dependent_couplings'] = \
                               self.write_print_parameters(self.coups_dep.values())


        file_py = self.read_template_file(self.param_template_py) % replace_dict
        
        return file_py

    def write_parameters(self, params, indent=8):
        """Write out the definitions of parameters"""

        # For each parameter type, write out the definition forcing a cast into the right type
        res_strings = []
        for param in params:
            # Not needed in Python, but kept for potential future use
            if param.type == 'real':
                res_strings.append('%sself.%s = %s(%s.real)'%(' '*indent, param.name, self.type_dict[param.type], param.name))                
            else:
                res_strings.append('%sself.%s = %s(%s)'%(' '*indent, param.name, self.type_dict[param.type], param.name))

        return '\n'.join(res_strings)

    def write_set_parameters(self, params, indent=8):
        """Write out the lines of independent parameters"""

        # For each parameter, write name = expr;

        res_strings = []
        for param in params:
            res_strings.append("%s%s" % (' '*indent, param.expr) )

        # Correct width sign for Majorana particles (where the width
        # and mass need to have the same sign)        
        for particle in self.model.get('particles'):
            if particle.is_fermion() and particle.get('self_antipart') and \
                   particle.get('width').lower() != 'zero':
                res_strings.append("%sif (%s < 0)" % (' '*indent,particle.get('mass')))
                res_strings.append("%(indent)s%(width)s = -abs(%(width)s)" % \
                                   {"width": particle.get('width'), 'indent': ' '*(indent+4)})

        return "\n".join(res_strings)

    def write_print_parameters(self, params, indent=8):
        """Write out the lines of independent parameters"""

        # For each parameter, write name = expr;

        res_strings = []
        for param in params:
            res_strings.append("%sres.append('{:<20s} = {:<20.16e}'.format('%s',self.%s))" % (' '*indent, param.name, param.name))
 
        return "\n".join(res_strings)

    # Routines for writing the ALOHA files

    def write_aloha_routines(self):
        """Generate the python aloha routines"""
       

        self.aloha_model.add_Lorentz_object(self.model.get('lorentz'))
        self.aloha_model.compute_subset(self.wanted_lorentz)
        # Write out the aloha routines in Python
        aloha_routines = []
        
        # First add the default external wavefunction routines
        wavefunction_routines = open(pjoin(MG5DIR,'aloha','template_files','wavefunctions.py'),'r').read() 
        open(pjoin(self.dir_path,'wavefunctions.py'),'w').write(
            'from __future__ import division\n'+wavefunction_routines)

        #aloha_routines.append(open(pjoin(MG5DIR,'aloha','template_files','wavefunctions.py'),'r').read())

        # Now write the process-depenent Feynman rules ones
        for routine in self.aloha_model.values():
            aloha_routines.append(routine.write(output_dir = None, 
                                                mode='mg5',
                                                language = 'Python'))

        for routine in self.aloha_model.external_routines:
            aloha_routines.append(
                     open(self.aloha_model.locate_external(routine, 'Python')).read())

        # Now collect imports
        python_imports = []
        new_aloha_routines = []
        for aloha_routine in aloha_routines:
            new_aloha_routine = []
            for line in aloha_routine.split('\n'):
                if any(line.startswith(token) for token in ['from','import']):
                    if line not in python_imports:
                        python_imports.append(line)
                else:
                    new_aloha_routine.append(line)
            new_aloha_routines.append('\n'.join(new_aloha_routine))
        aloha_routines = new_aloha_routines
       
        # Veto some imports
        vetoed_imports = ['import aloha.template_files.wavefunctions as wavefunctions']
        python_imports = [pi for pi in python_imports if pi not in vetoed_imports]
        python_imports.insert(0, 'import wavefunctions')

        aloha_output = open(pjoin(self.dir_path,'aloha_methods.py'),'w')
        aloha_output.write('from __future__ import division\n')
        # Write imports
        aloha_output.write('\n'.join(python_imports))
        aloha_output.write('\n'*2)

        # Write routines
        aloha_output.write('\n'.join(aloha_routines))

        aloha_output.close()

    #===============================================================================
    # Global helper methods
    #===============================================================================
    @classmethod
    def read_template_file(cls, filename, classpath=False):
        """Open a template file and return the contents."""
         
        return open(pjoin(plugin_path,'templates',filename),'r').read()

class ProcessOutputPython(export_v4.ProcessExporterFortranSA):
    check = False
    exporter = 'v4'
    output = 'dir'

class PythonMEExporter(export_python.ProcessExporterPython):

    matrix_method_template = 'matrix_method_python.py'

    def get_python_matrix_methods(self, gauge_check=False):
        """Write the matrix element calculation method for the processes"""

        replace_dict = {}

        # Extract version number and date from VERSION file
        info_lines = self.get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        for ime, matrix_element in enumerate(self.matrix_elements):
            process_string = matrix_element.get('processes')[0].shell_string()
            if process_string in self.matrix_methods:
                continue

            replace_dict['process_string'] = process_string

            # Extract number of external particles
            (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()
            replace_dict['nexternal'] = nexternal

            # Extract ncomb
            ncomb = matrix_element.get_helicity_combinations()
            replace_dict['ncomb'] = ncomb

            # Extract helicity lines
            helicity_lines = self.get_helicity_matrix(matrix_element)
            replace_dict['helicity_lines'] = helicity_lines

            # Extract overall denominator
            # Averaging initial state color, spin, and identical FS particles
            den_factor_line = self.get_den_factor_line(matrix_element)
            replace_dict['den_factor_line'] = den_factor_line

            # Extract process info lines for all processes
            process_lines = self.get_process_info_lines(matrix_element)
            replace_dict['process_lines'] = process_lines
        
            # Extract ngraphs
            ngraphs = matrix_element.get_number_of_amplitudes()
            replace_dict['ngraphs'] = ngraphs

            # Extract ndiags
            ndiags = len(matrix_element.get('diagrams'))
            replace_dict['ndiags'] = ndiags

            # Extract helas calls
            helas_calls = self.helas_call_writer.get_matrix_element_calls(\
                                                    matrix_element, gauge_check)
            replace_dict['helas_calls'] = "\n        ".join(helas_calls)

            # Extract nwavefuncs
            nwavefuncs = matrix_element.get_number_of_wavefunctions()
            replace_dict['nwavefuncs'] = nwavefuncs

            # Extract ncolor
            ncolor = max(1, len(matrix_element.get('color_basis')))
            replace_dict['ncolor'] = ncolor

            # Extract model parameter lines
            model_parameter_lines = \
                                 self.get_model_parameter_lines(matrix_element)
            replace_dict['model_parameters'] = model_parameter_lines

            # Extract color data lines
            color_matrix_lines = self.get_color_matrix_lines(matrix_element)
            replace_dict['color_matrix_lines'] = \
                                               "\n        ".join(color_matrix_lines)



            # External masses
            replace_dict['external_masses'] = self.get_external_masses(matrix_element)

            # Extract JAMP lines
            jamp_lines = self.get_jamp_lines(matrix_element)
            replace_dict['jamp_lines'] = "\n        ".join(jamp_lines)

            # Extract amp2 lines
            amp2_lines = self.get_amp2_lines(matrix_element,
                                        self.config_maps.setdefault(ime, []))
            replace_dict['amp2_lines'] = '\n        '.join(amp2_lines)

            method_file = open(os.path.join(plugin_path, 'templates',self.matrix_method_template)).read()
            method_file = method_file % replace_dict

            self.matrix_methods[process_string] = method_file

        return self.matrix_methods

    def get_external_masses(self, matrix_element):
        """ Return a tuple of the form (initial_masses, final_masses)."""

        
        proc_legs = matrix_element.get('processes')[0].get('legs')

        initial_masses = []
        final_masses = []
        for leg in proc_legs:
            if not leg.get('state'):
                initial_masses.append( self.model.get_particle(leg.get('id')).get('mass') )
            else:
                final_masses.append( self.model.get_particle(leg.get('id')).get('mass') )
        

        return '( (%s), (%s) )'%(
            ', '.join('model.%s'%mass for mass in initial_masses),
            ', '.join('model.%s'%mass for mass in final_masses)
        )

    def get_model_parameter_lines(self, matrix_element):
        """Return definitions for all model parameters used in this
        matrix element"""

        # Get all masses and widths used
        if aloha.complex_mass:
            parameters = [(wf.get('mass') == 'ZERO' or wf.get('width')=='ZERO') 
                          and wf.get('mass') or 'CMASS_%s' % wf.get('mass') 
                          for wf in \
                          matrix_element.get_all_wavefunctions()]
            parameters += [wf.get('mass') for wf in \
                      matrix_element.get_all_wavefunctions()]
        else:
            parameters = [wf.get('mass') for wf in \
                      matrix_element.get_all_wavefunctions()]
        parameters += [wf.get('width') for wf in \
                       matrix_element.get_all_wavefunctions()]
        parameters = list(set(parameters))
        if 'ZERO' in parameters:
            parameters.remove('ZERO')

        # Get all couplings used

        
        couplings = list(set([c.replace('-', '') for func \
                              in matrix_element.get_all_wavefunctions() + \
                              matrix_element.get_all_amplitudes() for c in func.get('coupling')
                              if func.get('mothers') ]))
        
        return "\n        ".join([\
                         "%(param)s = model.%(param)s"\
                         % {"param": param} for param in parameters]) + \
               "\n        " + "\n        ".join([\
                         "%(coup)s = model.%(coup)s"\
                              % {"coup": coup} for coup in couplings])

class PluginProcessExporterPython(object):
   
    exporter = 'v4'
    grouped_mode = False

    MEExporter=PythonMEExporter
    UFOModelConverter=UFOModelConverterPython
    check_sa_template = 'check_sa.py'

    def __init__(self, export_dir, helas_call_writers):
        self.export_dir = export_dir


        # Container to keep track of all_MEs exported
        self.all_MEs = []

        # Automatically add this output to the python path import system
        open(pjoin(self.export_dir,'__init__.py'),'w').write(
"""import sys
import os
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, root_path)
""")
        
        shutil.copytree( pjoin(plugin_path, 'templates','phase_space_generator'), 
                     pjoin(self.export_dir,'phase_space_generator') )
        os.makedirs(pjoin(self.export_dir, 'processes'))
        open(pjoin(self.export_dir,'processes','__init__.py'),'w')
        all_processes = open(pjoin(self.export_dir, 'processes','all_processes.py'),'w')
        # First add common imports
        all_processes.write('from __future__ import division\n')
        all_processes.write('from model.aloha_methods import *\n')
        all_processes.write('from model.wavefunctions import *\n')        
        all_processes.close()

        self.helas_call_writers = helas_call_writers

    def generate_subprocess_directory(self, matrix_element, dummy_helas_model, me_number):
        logger.info("Now generating Python output for %s"%(
                matrix_element.get('processes')[0].nice_string().replace('Process','process')))
        exporter = self.MEExporter(matrix_element, self.helas_call_writers)
        
        try:
            matrix_methods = exporter.get_python_matrix_methods(gauge_check=False)
            assert(len(matrix_methods)==1)
        except helas_call_writers.HelasWriterError, error:
            logger.critical(error)
            raise MadGraph5Error("Error when generation python matrix_element_methods.")

        self.all_MEs.extend(['Matrix_%s'%key for key in matrix_methods])

        all_processes = open(pjoin(self.export_dir, 'processes','all_processes.py'),'a')                            
        all_processes.write('\n'.join(matrix_methods.values()))
        all_processes.write('\n')
        all_processes.close()

        #for key, method in matrix_methods.items():
        #    open(pjoin(self.export_dir, 'processes','process_%s.py'%key),'w').write(method)

    def convert_model(self, model, wanted_lorentz, wanted_couplings):
        logger.info("Now outputting the model...")

        replace_dict = {}
        model_exporter = self.UFOModelConverter(
            model, pjoin(self.export_dir, 'model'), wanted_lorentz = wanted_lorentz,
            wanted_couplings = wanted_couplings, replace_dict=replace_dict)

        model_exporter.write_files()
        return

    def finalize(self, matrix_elements, history, options, flaglist):
        logger.info("Finalizing...")
        check_sa_path = pjoin(self.export_dir, 'check_sa.py')
        replace_dict = {'all_process_classes' : ', '.join(self.all_MEs)}
        open(check_sa_path,'w').write(
            open(pjoin(plugin_path, 'templates', self.check_sa_template),'r').read().format(**replace_dict))

        # Make check_sa.py executable
        os.chmod(check_sa_path, os.stat(check_sa_path).st_mode | stat.S_IEXEC)

# ==================================================================================================
# Classes below intended for the specialised TF output (which may or may not be needed in the end)
# ==================================================================================================

class ProcessOutputTF(ProcessOutputPython):
    check = False
    exporter = 'v4'
    output = 'dir'

class MEExporterTF(PythonMEExporter):
    matrix_method_template = 'matrix_method_TF.py'    
    # Specialisation for TF output to be implemented here
    pass

class UFOModelConverterTF(UFOModelConverterPython):
    param_template_py = 'model_template_TF.py'    
    # Specialisation for TF output to be implemented here
    pass 

class PluginProcessExporterTF(PluginProcessExporterPython):
    MEExporter=MEExporterTF
    UFOModelConverter=UFOModelConverterTF
    check_sa_template = 'check_sa_TF.py'    
    # Specialisation for TF output to be implemented here
    pass

class UFOHelasCallWriterTF(helas_call_writers.PythonUFOHelasCallWriter):
    # Specialisation for TF output to be implemented here    
    pass

