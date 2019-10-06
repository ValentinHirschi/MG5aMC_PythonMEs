#####################################################
#                                                   #
#  Source file of the PY8Kernels MG5aMC plugin.     #
#  Use only with consent of its author.             #
#                                                   #
#         author: Valentin Hirschi                  #
#                                                   #
#####################################################

import os
import logging
import itertools
import sys

from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
import madgraph.various.progressbar as pbar
import madgraph.interface.extended_cmd as cmd
import madgraph.iolibs.helas_call_writers as helas_call_writers
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc
import MG5aMC_PythonMEs.PluginExporters as PluginExporters

logger = logging.getLogger('MG5aMC_PythonMEs.Interface')

pjoin = os.path.join

class MG5aMC_PythonMEsPluginInterfaceError(MadGraph5Error):
    """ Error of the Exporter of the MG5aMC_PythonMEs interface. """

class MG5aMC_PythonMEsPluginInvalidCmd(InvalidCmd):
    """ Invalid command issued to the MG5aMC_PythonMEs interface. """

class MG5aMC_PythonMEsInterface(madgraph_interface.MadGraphCmd, cmd.CmdShell):
    """ Interface for steering the generation/output of MG5aMC_PythonMEs.
    We make it inherit from CmdShell so that launch_ext_prog does not attempt to start in WebMode."""
    
    def do_output(self, line):
        """ Wrapper to support the syntax output Python <args> """

        args = self.split_arg(line)
        if len(args)>=1 and args[0]=='Python':
            self.do_output_PythonMEs(' '.join(args[1:]))
        else:
            super(MG5aMC_PythonMEsInterface,self).do_output(' '.join(args))

    def do_output_PythonMEs(self, line):
        args = self.split_arg(line)
        super(MG5aMC_PythonMEsInterface,self).do_output(' '.join(['Python']+args))
        

    def export(self,*args,**opts):
        """Overwrite this so as to force a pythia8 type of output if the output mode is PY8MEs."""
        
        if self._export_format == 'plugin':
            # Also pass on the aloha model to the exporter (if it has been computed already)
            # so that it will be used when generating the model
            self._curr_exporter = PluginExporters.ProcessExporterPython(
                        self._export_dir,
                        helas_call_writers.PythonUFOHelasCallWriter(self._curr_model))

        super(MG5aMC_PythonMEsInterface,self).export(*args, **opts)

    #command to change the prompt 
    def preloop(self, *args, **opts):
        """only change the prompt after calling  the mother preloop command"""
        super(MG5aMC_PythonMEsInterface, self).preloop(*args,**opts)
        # The colored prompt screws up the terminal for some reason.
        #self.prompt = '\033[92mPY8Kernels > \033[0m'
        self.prompt = 'MG5aMC_PythonMEs > '
