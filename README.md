# TFMEs

This code is intended as a plugin to the High Energy Phsyics code [MadGraph5_aMC@NLO](https://launchpad.net/madgraph5).
It offers a new output mode for the standalone computation of Matrix Elements of scattering amplitudes. This special `TFMEs` format is in pure Python.

Over time, the goal is to adapt it so that it can be directly incorporated in a `TensorFlow` computational graph.

## Usage

Copy this project in the `PLUGIN` folder located in the root directory of your `MG5aMC` distribution (v2.6.6+).
The `MG5aMC` example script `test_TFMEs.mg5` can then simply be run as follows (from within the root directory of `MG5aMC`):
```
./bin/mg5_aMC --mode=TFMEs PLUGIN/TFMEs/test_TFMEs.mg5
```
The Python/TensorFlow code for this example selection of Matrix Elements will be generated in the folder `<MG5aMC_root_dir>/TFMEs_output_example` and its usage should be self-explanatory from reading the driver script at `<MG5aMC_root_dir>/TFMEs_output_example/check_sa.py`.
