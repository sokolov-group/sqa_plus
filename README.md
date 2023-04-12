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

### Matrix Block Evaluation: `sqa_plus.matrixBlock`
Automates the evaluation of terms applying normal-order to all active-space creation ahd annihilation operators with respect to physical vacuum, normal-order core creation and annihilation operators relative to the Fermi vaccum, and evaluate expectation values with respect to the active-space states.

```python
# Define spin-orbital Zeroth-order Heff
terms_Heff0 = sqa_plus.Heff(0)

# Evaluate terms
terms_Heff0 = matrixBlock(terms_Heff0)
```

### Export terms in NumPy's Einsum Notation: `sqa_plus.genEinsum`
Translates terms and SQA objects to [Numpy's Einsum](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html) notation, ready to be used in Python codes.
The `sqa_plus.genEinsum` function has options setted by positional arguments.

```python
sqa_plus.genEinsum(terms, lhs_str, ind_str, trans_rdm, trans_ind_str, suffix,
                          cvs_ind, val_ind, use_cvs_tensors, rm_trans_rdm_const, rm_core_int,
                          intermediate_list, opt_einsum_terms, optimize,
                          spin_integrated, use_spin_integrated_tensors,
                          help, **tensor_rename):
```

- `lhs_str`: Variable name to store the results of `np.einsum` operations;
- `ind_str`: String containing external indices of the summation;
- `suffix`: Add suffix to variables;
- `cvs_ind`: List of core indices to be included in Core-Valence Separation (CVS) space;
- `val_ind`: List of core indices not included in Core-Valence Separation (CVS) space;
- `use_cvs_tensors`: Generate expressions using CVS-ready variables;
- `rm_core_int`: Remove double-counted contributions to core terms when using the effective Hamiltonian;
- `opt_einsum_terms`: Select `optEinsum` or built-in Numpy `einsum` function;
- `optimize`: Configure `optimize` flag in `einsum`;
- `spin_integrated`: Express spin-integrated quantities by slices of spin-orbital objects;
- `use_spin_integrated_tensors`: Uses spin-integrated defined quantities instead, adding the spin cases as a suffix;


## Intermediates Generator
### Intermediates function: `sqa_plus.genIntermediates`
NOTE: Intermediates code requires [NumPy](https://numpy.org/) scientific computation package installed.

WARNING: This function has only been tested to function properly with the spin-orbital SQA terms. Work on integrating support for spin-adapted SQA terms is underway.

For a set of SQA terms, the function `sqa_plus.genIntermediates` is used to determine whether these terms contain expensive tensor contractions that are repeated multiple times. If so, the result of this contraction is defined and as an intermediate (`INT#`) tensor. Once all unique intermediate tensors are found, a new set of terms is defined with the intermediate tensors replacing the corresponding contraction in each original SQA term.  

```python
sqa_plus.genIntermediates(input_terms, ind_str = None, trans_rdm = False,
                          custom_path = None, factor_depth = 1, print_level = 4)
```
- `input_terms`: List of SQA term objects for which to generate intermediate tensors;
- `ind_str`: String containing external indices of the summation;
- `trans_rdm`: Boolean specifying whether either the bra or the ket of the expectation value is with respect to an excited state;
- `custom_path`: allows the user to specify (with zero-counting) a list of tuples that specify the order in which to contract the tensors in the terms. Not fully tested, use carefully;
- `factor_depth (default = 1)`: By default, the function will only generate intermediate tensors that are contractions of tensors found in the original set of SQA terms. If `factor_depth > 1`, the function will attempt to find new intermediates that it can define by re-using existing intermediate tensors;
- `print_level (default = 4)`: Control the level of detail that is being printed when searching through a set of terms to find intermediate tensors;

#### Example: Second-order contributions to M_00 block of the effective Hamiltonian matrix, M

$M^{(2)}_{ix,jy} = \langle \Psi_0 | [a^\dagger_i a_x, [\tilde{H}^{(2)}, a^\dagger_y a_j ]] | \Psi_0 \rangle $

```python
import sqa_plus as sqa

sqa.options.verbose = False # Default

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# External indices
i = sqa.index('I', [tg_c])
j = sqa.index('J', [tg_c])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])

# LHS
l_op  = [sqa.creOp(i), sqa.desOp(x)]
l_ind = 'IX'
# RHS
r_op  = [sqa.creOp(y), sqa.desOp(j)]
r_ind = 'JY'

# Define terms w/ desired operators
l_op_term = sqa.term(1.0, [], l_op)
r_op_term = sqa.term(1.0, [], r_op)

# Define Hamiltonian
print ("Creating effective Hamiltonian:\n")
effH = []
effH = sqa.Heff(order = 2, spin_integrated = False)

# First Commutator
print ("Performing first commutation...\n")
term1 = sqa.commutator(effH, r_op_term)

# Second Commutator
print ("Performing second commutation...\n")
term2 = sqa.commutator(l_op_term, term1)

# Evalauting terms
term3 = sqa.matrixBlock(term2)

# Filter through all terms keep only terms that involve 4-rdm
term4 = []
for term in term3:
    cre_des = 0

    for tens in term.tensors:
        if isinstance(tens, sqa.creOp) or isinstance(tens, sqa.desOp):
            cre_des += 1
        elif isinstance(tens, sqa.creDesTensor):
            cre_des += len(tens.ops)

    if cre_des == 8:
        term4.append(term)

# Original 4-RDM terms (no intermediates)
print ("Original 4-RDM terms, no intermediates:\n")
einsum_list_1 = sqa.genEinsum(term4, 'temp', l_ind + r_ind, rm_core_int = True)

print ("EINSUM CODE: 4-RDM terms, no intermediates:\n")
for einsum in einsum_list_1:
    print(einsum)
print('')
```

The 4-RDM terms for this matrix block are:

```
(  -0.25000) v(x,y,z,w) t1(I,x,u,v) t1(J,z,s,t) cre(Y) cre(y) cre(u) cre(v) des(X) des(w) des(s) des(t)
(  -0.25000) v(x,y,z,w) t1(I,x,u,v) t1(J,s,y,t) cre(Y) cre(u) cre(v) cre(s) des(X) des(z) des(w) des(t)  
(   0.25000) v(x,y,z,w) t1(I,x,u,v) t1(J,s,u,t) cre(Y) cre(y) cre(v) cre(s) des(X) des(z) des(w) des(t)
(   0.25000) v(x,y,z,w) t1(I,x,a,u) t1(J,v,a,s) cre(Y) cre(y) cre(u) cre(v) des(X) des(z) des(w) des(s)
(   0.25000) v(x,y,z,w) t1(I,u,x,v) t1(J,y,s,t) cre(Y) cre(z) cre(w) cre(v) des(X) des(u) des(s) des(t)
(  -0.12500) v(x,y,z,w) t1(I,u,x,v) t1(J,u,s,t) cre(Y) cre(z) cre(w) cre(v) des(X) des(y) des(s) des(t)
(   1.00000) v(x,y,z,w) t1(I,u,x,v) t1(J,s,z,t) cre(Y) cre(w) cre(v) cre(s) des(X) des(y) des(u) des(t)
(   0.25000) v(x,y,z,w) t1(I,u,x,v) t1(J,s,v,t) cre(Y) cre(z) cre(w) cre(s) des(X) des(y) des(u) des(t)
(  -0.25000) v(x,y,z,w) t1(I,u,a,x) t1(J,v,a,s) cre(Y) cre(z) cre(w) cre(v) des(X) des(y) des(u) des(s)
(  -0.12500) v(x,y,z,w) t1(I,i,x,u) t1(J,i,v,s) cre(Y) cre(z) cre(w) cre(u) des(X) des(y) des(v) des(s)
(   0.25000) v(x,y,z,w) t1(J,x,u,v) t1(I,s,u,t) cre(Y) cre(z) cre(w) cre(t) des(X) des(y) des(v) des(s)
(   0.25000) v(x,y,z,w) t1(J,x,a,u) t1(I,v,a,s) cre(Y) cre(z) cre(w) cre(s) des(X) des(y) des(u) des(v)
(  -0.12500) v(x,y,z,w) t1(J,u,x,v) t1(I,u,s,t) cre(Y) cre(y) cre(s) cre(t) des(X) des(z) des(w) des(v)
(   0.25000) v(x,y,z,w) t1(J,u,x,v) t1(I,s,v,t) cre(Y) cre(y) cre(u) cre(t) des(X) des(z) des(w) des(s)
(  -0.25000) v(x,y,z,w) t1(J,u,a,x) t1(I,v,a,s) cre(Y) cre(y) cre(u) cre(s) des(X) des(z) des(w) des(v)
(  -0.12500) v(x,y,z,w) t1(J,i,x,u) t1(I,i,v,s) cre(Y) cre(y) cre(v) cre(s) des(X) des(z) des(w) des(u)
(  -0.12500) kdelta(I,J) v(Y,x,y,z) t1(x,w,a,u) t1(v,s,a,t) cre(y) cre(z) cre(w) cre(t) des(X) des(u) des(v) des(s)
(   0.25000) kdelta(I,J) v(Y,x,y,z) t1(y,w,a,u) t1(v,s,a,t) cre(z) cre(u) cre(v) cre(s) des(X) des(x) des(w) des(t)
(   0.06250) kdelta(I,J) v(Y,x,y,z) t1(w,u,a,x) t1(v,s,a,t) cre(y) cre(z) cre(v) cre(s) des(X) des(w) des(u) des(t)
(  -0.12500) kdelta(I,J) v(Y,x,y,z) t1(w,u,a,y) t1(v,s,a,t) cre(z) cre(w) cre(u) cre(t) des(X) des(x) des(v) des(s)
(   0.06250) kdelta(I,J) v(Y,x,y,z) t1(i,x,w,u) t1(i,v,s,t) cre(y) cre(z) cre(s) cre(t) des(X) des(w) des(u) des(v)
(  -0.12500) kdelta(I,J) v(Y,x,y,z) t1(i,y,w,u) t1(i,v,s,t) cre(z) cre(w) cre(u) cre(v) des(X) des(x) des(s) des(t)
(  -0.12500) kdelta(I,J) v(Y,x,y,z) t1(i,w,x,u) t1(i,v,s,t) cre(y) cre(z) cre(u) cre(v) des(X) des(w) des(s) des(t)
(   0.25000) kdelta(I,J) v(Y,x,y,z) t1(i,w,y,u) t1(i,v,s,t) cre(z) cre(w) cre(s) cre(t) des(X) des(x) des(u) des(v)
(  -0.12500) kdelta(I,J) v(x,y,z,w) t1(Y,x,a,u) t1(v,s,a,t) cre(y) cre(u) cre(v) cre(s) des(X) des(z) des(w) des(t)
(   0.12500) kdelta(I,J) v(x,y,z,w) t1(Y,u,a,x) t1(v,s,a,t) cre(z) cre(w) cre(v) cre(s) des(X) des(y) des(u) des(t)
(  -0.12500) kdelta(I,J) v(x,y,z,w) t1(x,u,a,Y) t1(v,s,a,t) cre(z) cre(w) cre(u) cre(t) des(X) des(y) des(v) des(s)
(   0.25000) kdelta(I,J) v(x,y,z,w) t1(x,u,a,v) t1(Y,s,a,t) cre(z) cre(w) cre(u) cre(t) des(X) des(y) des(v) des(s)
(  -0.12500) kdelta(I,J) v(x,y,z,w) t1(x,u,a,v) t1(s,t,a,Y) cre(y) cre(v) cre(s) cre(t) des(X) des(z) des(w) des(u)
(  -0.12500) kdelta(I,J) v(x,y,z,w) t1(u,v,a,x) t1(Y,s,a,t) cre(y) cre(u) cre(v) cre(t) des(X) des(z) des(w) des(s)
(   0.06250) kdelta(I,J) v(x,y,z,w) t1(u,v,a,x) t1(s,t,a,Y) cre(z) cre(w) cre(s) cre(t) des(X) des(y) des(u) des(v)
(  -0.12500) kdelta(I,J) v(x,y,z,w) t1(i,Y,x,u) t1(i,v,s,t) cre(z) cre(w) cre(u) cre(v) des(X) des(y) des(s) des(t)
(   0.12500) kdelta(I,J) v(x,y,z,w) t1(i,x,Y,u) t1(i,v,s,t) cre(z) cre(w) cre(s) cre(t) des(X) des(y) des(u) des(v)
(   0.06250) kdelta(I,J) v(x,y,z,w) t1(i,x,u,v) t1(i,Y,s,t) cre(z) cre(w) cre(s) cre(t) des(X) des(y) des(u) des(v)
(  -0.12500) kdelta(I,J) v(x,y,z,w) t1(i,x,u,v) t1(i,s,Y,t) cre(y) cre(u) cre(v) cre(s) des(X) des(z) des(w) des(t)
(  -0.12500) kdelta(I,J) v(x,y,z,w) t1(i,u,Y,x) t1(i,v,s,t) cre(y) cre(u) cre(s) cre(t) des(X) des(z) des(w) des(v)
(  -0.12500) kdelta(I,J) v(x,y,z,w) t1(i,u,x,v) t1(i,Y,s,t) cre(y) cre(u) cre(s) cre(t) des(X) des(z) des(w) des(v)
(   0.25000) kdelta(I,J) v(x,y,z,w) t1(i,u,x,v) t1(i,s,Y,t) cre(z) cre(w) cre(v) cre(s) des(X) des(y) des(u) des(t)
```

The `sqa_plus.genEinsum` function can generate Python `np.einsum` calls to perform these contractions:  

```python
import np.einsum as einsum

temp =- 1/4 * einsum('xyzw,Ixuv,Jzst,YyuvXwst->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/4 * einsum('xyzw,Ixuv,Jsyt,YuvsXzwt->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('xyzw,Ixuv,Jsut,YyvsXzwt->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('xyzw,Ixau,Jvas,YyuvXzws->IXJY', v_aaaa, t1_caea, t1_caea, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('xyzw,Iuxv,Jyst,YzwvXust->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('xyzw,Iuxv,Just,YzwvXyst->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += einsum('xyzw,Iuxv,Jszt,YwvsXyut->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('xyzw,Iuxv,Jsvt,YzwsXyut->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/4 * einsum('xyzw,Iuax,Jvas,YzwvXyus->IXJY', v_aaaa, t1_caea, t1_caea, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('xyzw,Iixu,Jivs,YzwuXyvs->IXJY', v_aaaa, t1_ccaa, t1_ccaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('xyzw,Jxuv,Isut,YzwtXyvs->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('xyzw,Jxau,Ivas,YzwsXyuv->IXJY', v_aaaa, t1_caea, t1_caea, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('xyzw,Juxv,Iust,YystXzwv->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('xyzw,Juxv,Isvt,YyutXzws->IXJY', v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/4 * einsum('xyzw,Juax,Ivas,YyusXzwv->IXJY', v_aaaa, t1_caea, t1_caea, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('xyzw,Jixu,Iivs,YyvsXzwu->IXJY', v_aaaa, t1_ccaa, t1_ccaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,Yxyz,xwau,vsat,yzwtXuvs->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('IJ,Yxyz,ywau,vsat,zuvsXxwt->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/16 * einsum('IJ,Yxyz,wuax,vsat,yzvsXwut->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,Yxyz,wuay,vsat,zwutXxvs->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/16 * einsum('IJ,Yxyz,ixwu,ivst,yzstXwuv->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,Yxyz,iywu,ivst,zwuvXxst->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,Yxyz,iwxu,ivst,yzuvXwst->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('IJ,Yxyz,iwyu,ivst,zwstXxuv->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xyzw,Yxau,vsat,yuvsXzwt->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/8 * einsum('IJ,xyzw,Yuax,vsat,zwvsXyut->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xyzw,xuaY,vsat,zwutXyvs->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('IJ,xyzw,xuav,Ysat,zwutXyvs->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xyzw,xuav,staY,yvstXzwu->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xyzw,uvax,Ysat,yuvtXzws->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/16 * einsum('IJ,xyzw,uvax,staY,zwstXyuv->IXJY', np.identity(ncore), v_aaaa, t1_aaea, t1_aaea, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xyzw,iYxu,ivst,zwuvXyst->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/8 * einsum('IJ,xyzw,ixYu,ivst,zwstXyuv->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/16 * einsum('IJ,xyzw,ixuv,iYst,zwstXyuv->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xyzw,ixuv,isYt,yuvsXzwt->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xyzw,iuYx,ivst,yustXzwv->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xyzw,iuxv,iYst,yustXzwv->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
temp += 1/4 * einsum('IJ,xyzw,iuxv,isYt,zwvsXyut->IXJY', np.identity(ncore), v_aaaa, t1_caaa, t1_caaa, rdm_ccccaaaa, optimize = einsum_type)
```

Next, we can generate intermediates for these 4-RDM terms:

```python
# Scan for intermediates w/ factor_depth = 1
# INT terms will only be defined in terms of tensors in original equations
print ("Generating intermediate tensors:\n")
term5, intermediates = sqa.genIntermediates(term4, l_ind + r_ind, factor_depth = 1)

# Generate einsum for intermediates and new terms
int_einsum_list, einsum_list_2 = sqa.genEinsum(term5, 'temp', l_ind + r_ind, rm_core_int = True, intermediate_list = intermediates)

# Print intermediate term definitions
print ("EINSUM CODE: intermediate definitions:\n")
for int_def in int_einsum_list:
    print(int_def)
print('')

# Print terms w/ INT terms substitution
print ("EINSUM CODE: 4-RDM terms with INT tensor replacements:\n")
for einsum in einsum_list_2:
    print(einsum)
```

`sqa_plus.genIntermediates` makes a list of intermediate tensor definitions for all the unique and re-used contractions:

```python
INT1 = einsum('yYuvwXst,xyzw->YuvXstxz', rdm_ccccaaaa, v_aaaa, optimize = einsum_type)
INT2 = einsum('YuvszwXt,zwxy->YuvsXtxy', rdm_ccccaaaa, v_aaaa, optimize = einsum_type)
INT3 = einsum('yYvszwXt,xyzw->YvsXtx', rdm_ccccaaaa, v_aaaa, optimize = einsum_type)
INT4 = einsum('wYvsyXut,xyzw->YvsXutxz', rdm_ccccaaaa, v_aaaa, optimize = einsum_type)
INT5 = einsum('tyzwvsXu,vsat->yzwXua', rdm_ccccaaaa, t1_aaea, optimize = einsum_type)
INT6 = einsum('vuwXstzy,ivst->uwXzyi', rdm_ccccaaaa, t1_caaa, optimize = einsum_type)
```

After the list of intermediate tensor definitions, the original SQA terms are returned, with the INT tensors now in place of the tensor contractions that were initially in the term:

```python
temp =- 1/4 * einsum('Ixuv,Jzst,YuvXstxz->IXJY', t1_caaa, t1_caaa, INT1, optimize = einsum_type)
temp -= 1/4 * einsum('Ixuv,Jsyt,YuvsXtxy->IXJY', t1_caaa, t1_caaa, INT2, optimize = einsum_type)
temp -= 1/4 * einsum('Ixuv,Jsut,YvsXtx->IXJY', t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/4 * einsum('Ixau,Jvas,YuvXsx->IXJY', t1_caea, t1_caea, INT3, optimize = einsum_type)
temp -= 1/4 * einsum('Iuxv,Jyst,tsuXvYxy->IXJY', t1_caaa, t1_caaa, INT2, optimize = einsum_type)
temp += 1/8 * einsum('Iuxv,Just,tsXvYx->IXJY', t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp += einsum('Iuxv,Jszt,YvsXutxz->IXJY', t1_caaa, t1_caaa, INT4, optimize = einsum_type)
temp -= 1/4 * einsum('Iuxv,Jsvt,tuXsYx->IXJY', t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp += 1/4 * einsum('Iuax,Jvas,suXvYx->IXJY', t1_caea, t1_caea, INT3, optimize = einsum_type)
temp += 1/8 * einsum('Iixu,Jivs,svXuYx->IXJY', t1_ccaa, t1_ccaa, INT3, optimize = einsum_type)
temp -= 1/4 * einsum('Jxuv,Isut,svXtYx->IXJY', t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/4 * einsum('Jxau,Ivas,vuXsYx->IXJY', t1_caea, t1_caea, INT3, optimize = einsum_type)
temp += 1/8 * einsum('Juxv,Iust,YstXvx->IXJY', t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/4 * einsum('Juxv,Isvt,YutXsx->IXJY', t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp += 1/4 * einsum('Juax,Ivas,YusXvx->IXJY', t1_caea, t1_caea, INT3, optimize = einsum_type)
temp += 1/8 * einsum('Jixu,Iivs,YvsXux->IXJY', t1_ccaa, t1_ccaa, INT3, optimize = einsum_type)
temp += 1/8 * einsum('IJ,Yxyz,xwau,yzwXua->IXJY', np.identity(ncore), v_aaaa, t1_aaea, INT5, optimize = einsum_type)
temp -= 1/4 * einsum('IJ,Yxyz,ywau,wxXuza->IXJY', np.identity(ncore), v_aaaa, t1_aaea, INT5, optimize = einsum_type)
temp -= 1/16 * einsum('IJ,Yxyz,wuax,uwXzya->IXJY', np.identity(ncore), v_aaaa, t1_aaea, INT5, optimize = einsum_type)
temp += 1/8 * einsum('IJ,Yxyz,wuay,zwuXxa->IXJY', np.identity(ncore), v_aaaa, t1_aaea, INT5, optimize = einsum_type)
temp -= 1/16 * einsum('IJ,Yxyz,ixwu,uwXzyi->IXJY', np.identity(ncore), v_aaaa, t1_caaa, INT6, optimize = einsum_type)
temp += 1/8 * einsum('IJ,Yxyz,iywu,zwuXxi->IXJY', np.identity(ncore), v_aaaa, t1_caaa, INT6, optimize = einsum_type)
temp += 1/8 * einsum('IJ,Yxyz,iwxu,yzuXwi->IXJY', np.identity(ncore), v_aaaa, t1_caaa, INT6, optimize = einsum_type)
temp -= 1/4 * einsum('IJ,Yxyz,iwyu,uxXwzi->IXJY', np.identity(ncore), v_aaaa, t1_caaa, INT6, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,Yxau,vsat,uvsXtx->IXJY', np.identity(ncore), t1_aaea, t1_aaea, INT3, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,Yuax,vsat,tuXsvx->IXJY', np.identity(ncore), t1_aaea, t1_aaea, INT3, optimize = einsum_type)
temp += 1/8 * einsum('IJ,xuaY,vsat,svXtux->IXJY', np.identity(ncore), t1_aaea, t1_aaea, INT3, optimize = einsum_type)
temp -= 1/4 * einsum('IJ,xuav,Ysat,svXtux->IXJY', np.identity(ncore), t1_aaea, t1_aaea, INT3, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,xuav,staY,vstXux->IXJY', np.identity(ncore), t1_aaea, t1_aaea, INT3, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,uvax,Ysat,uvtXsx->IXJY', np.identity(ncore), t1_aaea, t1_aaea, INT3, optimize = einsum_type)
temp -= 1/16 * einsum('IJ,uvax,staY,vuXtsx->IXJY', np.identity(ncore), t1_aaea, t1_aaea, INT3, optimize = einsum_type)
temp += 1/8 * einsum('IJ,iYxu,ivst,tsXvux->IXJY', np.identity(ncore), t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,ixYu,ivst,vuXtsx->IXJY', np.identity(ncore), t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/16 * einsum('IJ,ixuv,iYst,vuXtsx->IXJY', np.identity(ncore), t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,ixuv,isYt,uvsXtx->IXJY', np.identity(ncore), t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,iuYx,ivst,ustXvx->IXJY', np.identity(ncore), t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/8 * einsum('IJ,iuxv,iYst,ustXvx->IXJY', np.identity(ncore), t1_caaa, t1_caaa, INT3, optimize = einsum_type)
temp -= 1/4 * einsum('IJ,iuxv,isYt,tuXsvx->IXJY', np.identity(ncore), t1_caaa, t1_caaa, INT3, optimize = einsum_type)
```

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
