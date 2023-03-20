# sqa_plus: SecondQuantizationAlgebra Plus

This program is intended to automate many of the tedious manipulations that one encounters when working in second quantization.
To use it, be sure that your ```$PYTHONPATH``` variable includes this directory.
Then, in the Python2 interpreter or in a Python2 script, use a command such as:

```python
import sqa_plus
```

# SecondQuantizationAlgebra Plus Capabilities

To see how to use core SecondQuantizationAlgebra capabilities, such as the definition of indexes, operators and terms, and algebraic functions, such as Normal Ordering, please visit the [original README file](sqa_original/README.md).

## Multireference Algebraic Diagrammatic Construction (MR-ADC) Theory Automation
### Spin-Orbital Effective Hamiltonian and Operators
Dyall Hamiltonian (`sqa_plus.dyallH`), Effective Hamiltonians up to second-order (`sqa_plus.Heff`), T amplitudes operators (`sqa_plus.getT`) and V perturbation operators (`sqa_plus.getV`) can be obtained automatically.

### Matrix Block Evaluation

## Intermediates Generator
### Intermediates function: `sqa_plus.genIntermediates`
Intermediates code requires [NumPy](https://numpy.org/) scientific computation package installed.

## Spin-Adapted Equation Generator
### Spin-Integrated Effective Hamiltonian and Operators
Spin-integrated effective Hamiltonians and operators can now be obtained by functions implemented in `sqaHeff.py`.
Dyall Hamiltonian (`sqa_plus.dyallH`), Effective Hamiltonians up to second-order (`sqa_plus.Heff`), T amplitudes operators (`sqa_plus.getT`) and V perturbation operators (`sqa_plus.getV`) are included.

```python
terms_Heff0_so = sqa_plus.Heff(0)                         # Spin-orbital zeroth-order Heff
terms_Heff0_si = sqa_plus.Heff(0, spin_integrated = True) # Spin-integrated zeroth-order Heff

terms_V_si = sqa_plus.getV(spin_integrated = True)
terms_T1_si = sqa_plus.getT(1, spin_integrated = True)
```

### Spin-Adaptation Automation from Spin-Integrated Equations: `sqa_plus.ConvertSpinIntegratedToSpinAdapted`
The function `sqa_plus.convertSpinIntegratedToAdapted` allows the automated conversion of spin-integrated terms to quantitites in the spin-adapted formulation.
Include the spin-adaptation of 1e- and 2e- integrals, and reduced density matrices (up to 4-RDMs).

#### Example
```python
# Define spin-integrated Zeroth-order Heff
terms_Heff0_si = sqa_plus.Heff(0, spin_integrated = True)

# Evaluate terms
terms_Heff0_si = matrixBlock(terms_Heff0_si)

# Convert to spin-adapted formulation
terms_Heff0_sa = sqa_plus.convertSpinIntegratedToAdapted(terms_Heff0_si)
```
The evaluated terms of the spin-integrated effective Hamiltonian in `terms_Heff0_si` are:
```
 (   1.00000) h(x,y) rdm(x,y)
 (   1.00000) h(x,y) rdm(x,y)
 (   1.00000) v(i,x,i,y) rdm(x,y)
 (   1.00000) v(i,x,i,y) rdm(x,y)
 (   1.00000) v(i,x,i,y) rdm(x,y)
 (   1.00000) v(i,x,i,y) rdm(x,y)
 (  -0.25000) v(x,y,z,w) rdm(x,y,w,z)
 (   0.25000) v(x,y,w,z) rdm(x,y,w,z)
 (  -0.25000) v(x,y,z,w) rdm(x,y,w,z)
 (  -0.25000) v(x,y,z,w) rdm(x,y,w,z)
 (  -0.25000) v(y,x,w,z) rdm(x,y,w,z)
 (   0.25000) v(y,x,z,w) rdm(x,y,w,z)
```

After the spin-adaptation procedure, `terms_Heff0_sa` should contain:
```
 (   1.00000) h(x,y) rdm(x,y)
 (   0.50000) v(x,y,z,w) rdm(x,y,z,w)
 (  -1.00000) v(x,i,i,y) rdm(x,y)
 (   2.00000) v(i,x,i,y) rdm(x,y)
```

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
