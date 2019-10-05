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
import TFMEs.TFMEsExporters as TFMEsExporters

logger = logging.getLogger('TFMEs_plugin.Interface')

pjoin = os.path.join

class TFMEsPluginInterfaceError(MadGraph5Error):
    """ Error of the Exporter of the TFMEs interface. """

class TFMEsPluginInvalidCmd(InvalidCmd):
    """ Invalid command issued to the TFMEs interface. """

class TFMEsInterface(madgraph_interface.MadGraphCmd, cmd.CmdShell):
    """ Interface for steering the generation/output of TFMEs.
    We make it inherit from CmdShell so that launch_ext_prog does not attempt to start in WebMode."""
    
    def do_output(self, line):
        """ Wrapper to support the syntax output TFMEs <args> """

        args = self.split_arg(line)
        if len(args)>=1 and args[0]=='TFMEs':
            self.do_output_TFMEs(' '.join(args[1:]))
        else:
            super(TFMEsInterface,self).do_output(' '.join(args))

    def do_output_TFMEs(self, line):
        args = self.split_arg(line)
        super(TFMEsInterface,self).do_output(' '.join(['TFMEs']+args))
        

    def export(self,*args,**opts):
        """Overwrite this so as to force a pythia8 type of output if the output mode is PY8MEs."""
        
        if self._export_format == 'plugin':
            # Also pass on the aloha model to the exporter (if it has been computed already)
            # so that it will be used when generating the model
            self._curr_exporter = TFMEsExporters.ProcessExporterTF(
                        self._export_dir,
                        helas_call_writers.PythonUFOHelasCallWriter(self._curr_model))

        super(TFMEsInterface,self).export(*args, **opts)

    #command to change the prompt 
    def preloop(self, *args, **opts):
        """only change the prompt after calling  the mother preloop command"""
        super(TFMEsInterface, self).preloop(*args,**opts)
        # The colored prompt screws up the terminal for some reason.
        #self.prompt = '\033[92mPY8Kernels > \033[0m'
        self.prompt = 'TFMEs > '
