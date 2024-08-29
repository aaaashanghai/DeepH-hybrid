# DeepH-hybrid

This code provides an implementation of the DeepH-hybrid method described in the paper _A deep equivariant neural network approach for efficient hybrid density functional calculations_ ([arXiv:2302.08221](https://arxiv.org/abs/2302.08221)).

DeepH-hybrid is built upon DeepH series of Deep-learning Hamiltonians. The case studies in DeepH-hybrid is carried out with [DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3), but may also be integrated with [DeepH-pack](https://github.com/mzjb/DeepH-pack). The current code contains only the newly included code of DeepH-hybrid, including a brief demo for these parts.

## Installation
The newly added code is executable under the same python environment of DeepH-E3. In case of running this code under a different environment, a Python > 3.9 is required, with following packages installed:

- numpy
- h5py

## Usage
The usage of DeepH-hybrid is in full analogy with [DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3). However, here's a few things to note:

### Preparation of DFT-hybrid Data
Preparation of DFT-hybrid data for training DeepH-hybrid involves using the [ABACUS](http://abacus.deepmodeling.com/en/latest/) package. The [ABACUS interface for DeepH](https://github.com/mzjb/DeepH-pack/blob/main/deeph/preprocess/abacus_get_data.py), initially developed for DeepH-hybrid, is already open-sourced in the [DeepH-pack](https://github.com/mzjb/DeepH-pack).

For usage of ABACUS, please refer to [the manual of the ABACUS package](http://abacus.deepmodeling.com/en/latest/). For interfacing ABACUS with DeepH-hybrid, the `out_mat_hs2` option of the ABACUS must be set to `1`.

### Modifying the $R_{\text{C}}^{\text{hyb}}$
The modification of the cutoff radius, $R_{\text{C}}^{\text{hyb}}$, is done before the training process. A config file must be specified for determining $R_{\text{C}}^{\text{hyb}}$, with an example in `example/input/rc_config.json`. The config is stored in a json file. Each key and value of the json files specifies the atomic indices (1 for H, 2 for He, etc.) and the cutoff radius for this element (length unit in Bohr). The demo takes the value of $R_{\text{C}}^{\text{hyb}}=2R_{\text C}$, which is used in cased studies of the main text.

A simple example use is provided in the `example` directory. Just execute `bash run.sh` at that directory, and the `example/input` directory will be modified into files in the `example/output` directory.

The current version of [DeepH-E3](https://github.com/Xiaoxun-Gong/DeepH-E3) supports DeepH-hybrid usage after this modification is accomplished between `deeph-preprocess` and `deeph-train`.