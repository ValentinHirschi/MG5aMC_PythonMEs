# MG5aMC_PythonMEs

This code is intended as a plugin to the High Energy Phsyics code [MadGraph5_aMC@NLO](https://launchpad.net/madgraph5) (v2.6.6+).
It offers a new output mode for the standalone computation of Matrix Elements of scattering amplitudes. This special `MG5aMC_PythonMEs` format is in pure Python.

## Usage

Copy this project in the `PLUGIN` folder located in the root directory of your `MG5aMC` distribution (v2.6.6+).
The `MG5aMC` example script `test_MG5aMC_PythonMEs.mg5` can then simply be run as follows (from within the root directory of `MG5aMC`):
```
./bin/mg5_aMC --mode=MG5aMC_PythonMEs PLUGIN/MG5aMC_PythonMEs/test_MG5aMC_PythonMEs.mg5
```
The Python/TensorFlow code for this example selection of Matrix Elements will be generated in the folder `<MG5aMC_root_dir>/MG5aMC_PythonMEs_output_example` and its usage should be self-explanatory from reading the driver script at `<MG5aMC_root_dir>/MG5aMC_PythonMEs_output_example/check_sa.py`.

## Further Development

Over time, the goal is to adapt this plugin so that its output can be directly incorporated in ML engines, such as in a `TensorFlow` computational graph.

To this effect the current repository provides the additional format `TF` for the `output` command.
For now, this output behaves exactly identically to the `Python` output format (see by yourself by running `test_MG5aMC_TFMEs.mg5`), but uses separate daughter classes and templates (suffixed with `TF`) which are ready to be specialised as needed for the output code to be directly inserted in `TensorFlow`. In particular the driving script template `check_sa_TF.py` should be turned into a demo `Jupyter` notebook demonstrating the integration within `TensorFlow`.

Note that for now the model independent parameters are hard-coded to their default value in the `parameters.py` script. Eventually we may want to add a facility for reading an SLHA input card, but this is not needed for now.
