The .py files in this directory are those used to generate the
terms involved in the various commutator and residual functions
involved in Canonical Transformation theory.  Before running or
modifying these scripts, we suggest looking at the tutorial.

The .out files in this directory contain the output for the
various scripts.  The final list of resulting terms is at the
end of the file.  Below is a description of the tensor names
used in the scripts and output files.

h1aa     -- all-alpha spin 1-body integrals
h1bb     -- all-beta spin 1-body integrals
h2aaaa   -- all-alpha spin 2-body integrals
h2bbbb   -- all-beta spin 2-body integrals
h2baba   -- mixed spin 2-body integrals
a1aa     -- all-alpha spin 1-body anti-symmetric amplitudes
a1bb     -- all-beta spin 1-body anti-symmetric amplitudes
a2aaaa   -- all-alpha spin 2-body anti-symmetric amplitudes
a2bbbb   -- all-beta spin 2-body anti-symmetric amplitudes
a2baba   -- mixed spin 2-body anti-symmetric amplitudes
d1aa     -- all-alpha spin 1-body reduced density matrix
d1bb     -- all-beta spin 1-body reduced density matrix
d2aaaa   -- all-alpha spin 2-body reduced density matrix
d2bbbb   -- all-beta spin 2-body reduced density matrix
d2baba   -- mixed spin 2-body reduced density matrix
d3aaaaaa -- all-alpha spin 3-body reduced density matrix
d3bbbbbb -- all-beta spin 3-body reduced density matrix
d3abbabb -- mixed spin 3-body reduced density matrix
d3baabaa -- mixed spin 3-body reduced density matrix

For commutator functions, r is used to label the transformed
integrals for the new effective hamiltonian term H_bar = [H, A].
Each term has an r tensor, which is used to store the result of
the contractions of the other tensors.

r0       -- constant term for the transformed integrals
r1aa     -- all-alpha spin 1-body transformed integrals
r1bb     -- all-beta spin 1-body transformed integrals
r2aaaa   -- all-alpha spin 2-body transformed integrals
r2bbbb   -- all-beta spin 2-body transformed integrals
r2baba   -- mixed spin 2-body transformed integrals

For residual functions, r is used to label the 1- and 2-body
residuals.  Again, the r tensor is used to store the result of
the contractions of the other tensors.

r1aa     -- all-alpha spin 1-body residuals
r1bb     -- all-beta spin 1-body residuals
r2aaaa   -- all-alpha spin 2-body residuals
r2bbbb   -- all-beta spin 2-body residuals
r2baba   -- mixed spin 2-body residuals
