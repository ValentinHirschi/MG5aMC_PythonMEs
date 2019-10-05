## import the required files
import os
import sys
root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.insert(0, root_path)

import TFMEs.TFMEsInterface as TFMEsInterface
import TFMEs.TFMEsExporters as TFMEsExporters

# Three types of functionality are allowed in a plugin
#   1. new output mode
#   2. new cluster support
#   3. new interface

# 1. Define new output mode
#    example: new_output = {'myformat': MYCLASS}
#    madgraph will then allow the command "output myformat PATH"
#    MYCLASS should inherated of the class madgraph.iolibs.export_v4.VirtualExporter 
new_output = {'TFMEs': TFMEsExporters.ProcessOutputTF}

# 2. Define new way to handle the cluster.
#    example new_cluster = {'mycluster': MYCLUSTERCLASS}
#    allow "set cluster_type mycluster" in madgraph
#    MYCLUSTERCLASS should inherated from madgraph.various.cluster.Cluster
new_cluster = {}

# 3. Define a new interface (allow to add/modify MG5 command)
#    This can be activated via ./bin/mg5_aMC --mode=PLUGINNAME
## Put None if no dedicated command are required
new_interface = TFMEsInterface.TFMEsInterface
 
 
########################## CONTROL VARIABLE ####################################
__author__ = 'Valentin Hirschi'
__email__ = 'valentin.hirschi@gmail.com'
__version__ = (1,0,0)
minimal_mg5amcnlo_version = (2,6,6) 
maximal_mg5amcnlo_version = (1000,1000,1000)
latest_validated_version = (2,6,6)
