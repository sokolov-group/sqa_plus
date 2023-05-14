# sqa_plus: SecondQuantizationAlgebra Plus

This program is intended to automate many of the tedious manipulations that one encounters when working in second quantization.
To use it, be sure that your ```$PYTHONPATH``` variable includes this directory.
Then, in the Python2 interpreter or in a Python2 script, use a command such as:

```python
import sqa_plus
```

# SecondQuantizationAlgebra Plus Capabilities

To see how to use core SecondQuantizationAlgebra capabilities, such as the definition of indexes, operators and terms, and algebraic functions, such as Normal Ordering, please visit the [original README file](README_sqa_extra.md).

# Multireference Algebraic Diagrammatic Construction (MR-ADC) Theory Automation
## Spin-Orbital Effective Hamiltonian and Operators
Dyall Hamiltonian (`sqa_plus.dyallH`), Effective Hamiltonians up to second-order (`sqa_plus.Heff`), T amplitudes operators (`sqa_plus.Tamplitude`) and V perturbation operators (`sqa_plus.Vperturbation`) can be obtained automatically.

## Spin-Integrated Effective Hamiltonian and Operators
Spin-integrated effective Hamiltonians and operators can now be obtained by functions implemented in `sqaHeff.py`.
Dyall Hamiltonian (`sqa_plus.dyallH`), Effective Hamiltonians up to second-order (`sqa_plus.Heff`), T amplitudes operators (`sqa_plus.Tamplitude`) and V perturbation operators (`sqa_plus.Vperturbation`) are included.
The spin-integrated basis can be globally selected by using `sqa_plus.options.spin_integrated = True`.
By default, `sqa_plus` assumes the spin-orbital basis (`sqa_plus.options.spin_orbital = True`).

```python
sqa_plus.options.spin_integrated = True
terms_Heff0_si = sqa_plus.Heff(0)        # Spin-integrated zeroth-order Heff
```

## Core-Valence Separation (CVS) approximattion using SQA+ Equation Generator
In a similar design to the choice of the spin basis, the CVS approximation can be selected by using `sqa_plus.options.cvs_approach = True`.
It will set globally all the implemented operators (`sqaHeff.py`, `sqa_plus.dyallH`, `sqa_plus.Heff`, `sqa_plus.Tamplitude`, and `sqa_plus.Vperturbation`) to assume the CVS approximation.

In the SQA+ language, the CVS approximation decouple the `options.core_type` type indices in two new types of indices: `options.cvs_core_type` and `options.cvs_valence_type`.
Then, core indices defined by the user must be of these types.

## Using summed (dummy) indices lists: `sqa_plus.indexLists`
Besides the usage of `sqa_plus.index` class to define dummy indices, we provide an automated way to generate unique indices.
The `sqa_plus.indexLists` object can be used to create dummy indices, free of conflict with SQA+ implemented operators.
New indices can be obtained using the `new_index()` method.

```python
# Examples of sqa_plus.indexLists usage
cvs_core_ind_1 = sqa_plus.indexLists.cvs_core.new_index()
active_ind_1 = sqa_plus.indexLists.active.new_index()
virtual_alpha_ind_1 = sqa_plus.indexLists.virtual_alpha.new_index()
```

#### Dummy index types supported
- Spin-orbital basis:
	- `indexLists.core`
	- `indexLists.active`
	- `indexLists.virtual`

- Spin-orbital basis using CVS:
	- `indexLists.cvs_core`
	- `indexLists.cvs_valence`
	- `indexLists.active`
	- `indexLists.virtual`

- Spin-integrated basis:
	- `indexLists.core_alpha`
	- `indexLists.active_alpha`
	- `indexLists.virtual_alpha`

	- `indexLists.core_beta`
	- `indexLists.active_beta`
	- `indexLists.virtual_beta`

- Spin-integrated basis using CVS:
	- `indexLists.cvs_core_alpha`
	- `indexLists.cvs_valence_alpha`
	- `indexLists.active_alpha`
	- `indexLists.virtual_alpha`

	- `indexLists.cvs_core_beta`
	- `indexLists.cvs_valence_beta`
	- `indexLists.active_beta`
	- `indexLists.virtual_beta`

## Spin-Adapted Equation Generator
### Spin-Adaptation Automation from Spin-Integrated Equations: `sqa_plus.ConvertSpinIntegratedToSpinAdapted`
The function `sqa_plus.convertSpinIntegratedToAdapted` allows the automated conversion of spin-integrated terms to quantitites in the spin-adapted formulation.
Include the spin-adaptation of 1e- and 2e- integrals, and reduced density matrices (up to 4-RDMs).

#### Example
```python
# Define spin-integrated Zeroth-order Heff
sqa_plus.options.spin_integrated = True
terms_Heff0_si = sqa_plus.Heff(0)

# Evaluate terms
terms_Heff0_si = sqa_plus.matrixBlock(terms_Heff0_si)

# Convert to spin-adapted formulation
terms_Heff0_sa = sqa_plus.convertSpinIntegratedToAdapted(terms_Heff0_si)
```
The evaluated terms of the spin-integrated effective Hamiltonian in `terms_Heff0_si` are:
```
 (   1.00000) E_fc(Const.) 
 (   1.00000) h(x,y) cre(x) des(y) 
 (   1.00000) h(x,y) cre(x) des(y) 
 (   1.00000) v(i,x,i,y) cre(x) des(y) 
 (   1.00000) v(i,x,i,y) cre(x) des(y) 
 (   1.00000) v(i,x,i,y) cre(x) des(y) 
 (   1.00000) v(i,x,i,y) cre(x) des(y) 
 (  -0.25000) v(x,y,z,w) cre(x) cre(y) des(z) des(w) 
 (  -1.00000) v(x,y,z,w) cre(x) cre(y) des(z) des(w) 
 (  -0.25000) v(x,y,z,w) cre(x) cre(y) des(z) des(w) 
```

After the spin-adaptation procedure, `terms_Heff0_sa` should contain:
```
 (   1.00000) E_fc(Const.) 
 (   1.00000) h(x,y) rdm(x,y) 
 (   2.00000) v(i,i,x,y) rdm(x,y) 
 (  -1.00000) v(i,x,y,i) rdm(x,y) 
 (   0.50000) v(x,y,z,w) rdm(x,z,y,w) 
```

Note that `sqa_plus.convertSpinIntegratedToAdapted` will globally set `sqa_plus.options.spin_adapted = True`) and changes 2e- integrals and RDMs to the chemists' notation (`sqa_plus.options.chemists_notation = True`).

## Matrix Block Evaluation: `sqa_plus.matrixBlock`
Automates the evaluation of terms applying normal-order to all active-space creation ahd annihilation operators with respect to physical vacuum, normal-order core creation and annihilation operators relative to the Fermi vaccum, and evaluate expectation values with respect to the active-space states.

```python
# Define spin-orbital Zeroth-order Heff
terms_Heff0 = sqa_plus.Heff(0)

# Evaluate terms
terms_Heff0 = matrixBlock(terms_Heff0)
```

## Intermediates Generator
### Intermediates function: `sqa_plus.genIntermediates`
Intermediates code requires [NumPy](https://numpy.org/) scientific computation package installed.

## Export terms in NumPy's Einsum Notation: `sqa_plus.genEinsum`
Translates terms and SQA+ objects to [Numpy's Einsum](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html) notation, ready to be used in Python codes.
The `sqa_plus.genEinsum` function has options setted by positional arguments.

```python
sqa_plus.genEinsum(terms, lhs_string = None, indices_string = None, suffix = None,
				   trans_indices_string = None, intermediate_list = None, help = False, **tensor_rename):
```

- `lhs_string`: Variable name to store the results of `np.einsum` operations;
- `indices_string`: String containing external indices of the summation;
- `suffix`: Add suffix to variables;

Aditionally, options for genEinsum can be selected using the object `sqa_plus.options`

- `sqa_plus.options.genEinsum.lhs_string`: Variable name to store the results of `np.einsum` operations;
- `sqa_plus.options.genEinsum.indices_string`: String containing external indices of the summation;

- `sqa_plus.options.genEinsum.spin_integrated_tensors`: Express spin-integrated quantities by slices of spin-orbital objects; (default: `False`)
- `sqa_plus.options.genEinsum.keep_user_defined_dummy_names`: (default: `False`)

- `sqa_plus.options.genEinsum.trans_rdm`: Use transition reduced density matrices (tRDMs); (default: `False`)
- `sqa_plus.options.genEinsum.trans_indices_string`: String containing tRDMs indices;
- `sqa_plus.options.genEinsum.suffix`: Add suffix to variables;

- `sqa_plus.options.genEinsum.cvs_indices_list`: List of core indices to be included in Core-Valence Separation (CVS) space;
- `sqa_plus.options.genEinsum.valence_indices_list`: List of core indices not included in Core-Valence Separation (CVS) space;
- `sqa_plus.options.genEinsum.cvs_tensors`: Generate expressions using CVS-ready variables; (default: `True`)

- `sqa_plus.options.genEinsum.remove_trans_rdm_constant`: Remove tRDMs constants (default: `False`)
- `sqa_plus.options.genEinsum.remove_core_integrals`: Remove double-counted contributions to core terms when using the effective Hamiltonian (default: `True`)
- `sqa_plus.options.genEinsum.intermediate_list`: List of intermediates tensors;

- `sqa_plus.options.genEinsum.opt_einsum_terms`: Select `optEinsum` or built-in Numpy `einsum` function; (default: `True`)
- `sqa_plus.options.genEinsum.optimize`: Configure `optimize` flag in `einsum`;(default: `True`)

# Credits
## SecondQuantizationAlgebra
SecondQuantizationAlgebra was originally developed by Eric Neuscamman <eric.neuscamman@gmail.com>

In addition, any modification or use of this software should cite the following paper:
```
  E. Neuscamman, T. Yanai, and G. K.-L. Chan.
  J. Chem. Phys. 130, 124102 (2009)
```

## SecondQuantizationAlgebra Plus
### MR-ADC Automation
- Alexander Yu. Sokolov <alexander.y.sokolov@gmail.com>
- Koushik Chatterjee <koushikchatterjee7@gmail.com>

### Intermediates Generator
- Ilia Mazin <ilia.mazin@gmail.com>

### Spin-Adaptation Automation
- Carlos E. V. de Moura <carlosevmoura@gmail.com>

### Core-Valence Separation Approximation
- Carlos E. V. de Moura <carlosevmoura@gmail.com>
- Ilia Mazin <ilia.mazin@gmail.com>
