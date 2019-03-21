# file:  sqaAddon.py
# author:  Koushik Chatterjee
# date:  August 31, 2018
#
# summary:
#           matrixBlock: General routine used to construct matrix block.
#           dummyLabel: Relabel dummy indices according to core, active, external.
#           filterVirtual: Calculate expectation value wrt virtual orbitals.
#           filterCore: Calculate expectation value wrt core orbitals.
#           normalOrderCore: Normal order with respect to core orbitals.
#
# Copyright (C) 2018-2019 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
#

from sqaIndex import index
from sqaTensor import tensor, kroneckerDelta, sfExOp, creOp, desOp, creDesTensor
from sqaTerm import term, multiplyTerms, termChop, sortOps, combineTerms
from sqaOptions import options
from sqaMisc import makeTuples, allDifferent, makePermutations
from sqaSymmetry import symmetry

from sqaNormalOrder import normalOrder

# Create an object of the index class
#addon = index()

#####################################
def matrixBlock(nterms, transRDM = False):
 "Construct matrix block."
#
 print ""
 print "################ Addon ################"
 print_header()
 print ""
 fTerms = []
# Dummy indices label upate
# dummyLabel(nterms)
#
# Filter zero terms wrt virtual (note: Filter first for virtual orbitals)
 filterVirtual(nterms)
 termChop(nterms)
#
# Normal ordering with respect to core orbitals
# fTerms = []
 for t in nterms:
#    print 'Term=', t
    trm = normalOrderCore(t)
    fTerms.extend(trm)
#
# Evaluate Kroneker delta
 for t in fTerms:
    t.contractDeltaFuncs()
#
# termChop(fTerms)
#
# Filter zero terms wrt core (note: after filtering virtual, then do for core)
 filterCore(fTerms)
 termChop(fTerms)
#
 combineTerms(fTerms)
#
## Dummy indices label upate
# dummyLabel(fTerms)
#
# Contract delta function for both non-dummy indices
 contractDeltaFuncs_nondummy(fTerms)
#
# If (transRDM = True) =>  Remove those constant terms
 if (transRDM):
    for trm in fTerms:
        iremove = True
        for i in range(len(trm.tensors)):
            t = trm.tensors[i]
            if (isinstance(t, creOp) or isinstance(t, desOp)):
               iremove = False
        if (iremove):
            trm.numConstant = 0.0
    termChop(fTerms)
#
# Dummy indices label upate
 dummyLabel(fTerms)
#
# Reorder tensor indices: (core < active < virtual) order
 tensorIndex_order(fTerms)
#
# Print the final results
 print ""
 print "####### Final results:#######"
 for t in fTerms:
    index_types = ()
    for t_tensor in t.tensors:
        for t_tensor_index in t_tensor.indices:
            index_types += t_tensor_index.indType[0]
    print t
#    print index_types
#
 print ""
 return fTerms
#
#####################################
def dummyLabel(nterms):
 "A function to relabel dummy indices."
#
 print "Dummy indices label update:=>"
##
 for t in nterms:
    mymap ={}
    index_types = ()
#
    coreInd = list('ijklmnopq')
    actvInd = list('xyzwuvstr')
    virtInd = list('abcdefgh')
#
    for t_tensor in t.tensors:
#
        for t_tensor_index in range(len(t_tensor.indices)):
#
            # Decide which new label to assign
            index_type = t_tensor.indices[t_tensor_index].indType
            index_name = t_tensor.indices[t_tensor_index].name
            index_summed = t_tensor.indices[t_tensor_index].isSummed
            if index_summed:
                if len(index_type) > 1:
                    raise Exception('This code does not support general indices for now')
                else:
                    if index_name not in mymap.keys():
                        if (index_type[0][0] == 'core'):
                            mymap[index_name] = coreInd[0]
                            coreInd.pop(0)
                        elif (index_type[0][0] == 'active'):
                            mymap[index_name] = actvInd[0]
                            actvInd.pop(0)
                        else:
                            mymap[index_name] = virtInd[0]
                            virtInd.pop(0)

                    # Update the label
                    t_tensor.indices[t_tensor_index].name = mymap[index_name]
#
#        index_types += t_tensor.indices[t_tensor_index].indType
    print t
#    print t.numConstant
#    print index_types
#
#####################################
#
def filterVirtual(nterms):
 "A function to calculate expectation value wrt virtual: filter zero terms wrt virtual."
#
 print ""
 print "Expectation value: Filter zero terms wrt virtual:=>"
##
 for t in nterms:
    mymap ={}
    index_types = ()
#
    for t_tensor in t.tensors:
    #     print t_tensor, t_tensor.name, t_tensor.indices
         for t_tensor_index in range(len(t_tensor.indices)):
#
              index_type = t_tensor.indices[t_tensor_index].indType
              index_name = t_tensor.indices[t_tensor_index].name
              tensor_name = t_tensor.name
#
              if len(index_type) > 1:
                 raise Exception('This code does not support general indices for now')
              else:
#
                 if (tensor_name == 'des'):
                      if (index_type[0][0] == 'virtual'):
                         t.numConstant = 0.0
#                         print 'cond1:',t_tensor, t_tensor.name, t_tensor.indices[t_tensor_index].name,t_tensor.indices[t_tensor_index].indType[0][0]
                 elif (tensor_name == 'cre'):
                      if (index_type[0][0] == 'virtual'):
                        t.numConstant = 0.0
#                         print 'cond2:',t_tensor, t_tensor.name, t_tensor.indices[t_tensor_index].name,t_tensor.indices[t_tensor_index].indType[0][0]
#
#
    print t
#
#####################################
#
def filterCore(nterms):
 "A function to calculate expectation value wrt core: filter zero terms wrt core."
#
 print ""
 print "Expectation value: Filter zero terms wrt core:=>"
##
 for t in nterms:
    mymap ={}
    index_types = ()
#
    for t_tensor in t.tensors:
    #     print t_tensor, t_tensor.name, t_tensor.indices
         for t_tensor_index in range(len(t_tensor.indices)):
#
              index_type = t_tensor.indices[t_tensor_index].indType
              index_name = t_tensor.indices[t_tensor_index].name
              tensor_name = t_tensor.name
#
              if len(index_type) > 1:
                 raise Exception('This code does not support general indices for now')
              else:
#
                 if (tensor_name == 'des'):
                      if (index_type[0][0] == 'core'):
                         t.numConstant = 0.0
#                         print 'cond1:',t_tensor, t_tensor.name, t_tensor.indices[t_tensor_index].name,t_tensor.indices[t_tensor_index].indType[0][0]
                 elif (tensor_name == 'cre'):
                      if (index_type[0][0] == 'core'):
                        t.numConstant = 0.0
#                         print 'cond2:',t_tensor, t_tensor.name, t_tensor.indices[t_tensor_index].name,t_tensor.indices[t_tensor_index].indType[0][0]
#
#
    print t
#
#####################################
#
def normalOrderCore(inTerm):
  "Returns a list of terms resulting from normal ordering the operators in inTerm."
  print ""
  print "Normal ordering with respect to core:=>"
  print 'Term=', inTerm
#  if options.verbose:
#    print "converting to normal order:  %s" %(str(inTerm))

  # check that inTerm is a term
  if not isinstance(inTerm, term):
    raise TypeError, "inTerm must be of class term"

  # determine what types of operators the term contains
  has_creDesOps = False
  has_sfExOps = False
  for t in inTerm.tensors:
    if isinstance(t, creOp) or isinstance(t, desOp):
      has_creDesOps = True
    elif isinstance(t, sfExOp):
      has_sfExOps = True

  # If term has both creation/destruction operators and spin free excitation operators,
  # raise an error
  if has_creDesOps and has_sfExOps:
    raise RuntimeError, "Normal ordering not implemented when both creOp/desOp and sfExOp " + \
                        "tensors are present"

  # if the term is already normal ordered, return it unchanged
 # elif inTerm.isNormalOrdered():
 #   print "koushik ordering"
 #   outTerms = [inTerm.copy()]

  # Normal ordering for creOp/desOp
  elif has_creDesOps:

    # Separate the cre/des operators from other tensors
    ops = []
    nonOps = []
    for t in inTerm.tensors:
      if isinstance(t, creOp) or isinstance(t, desOp):
        ops.append(t.copy())
      else:
        nonOps.append(t.copy())

    # Generate all contraction pairs
##koushik
    contractionPairs = []
    for i in range(len(ops)):
       iTerm = ops[i]
       iType = iTerm.indices[0].indType
       for j in range(i+1,len(ops)):
#         if (jTerm.indices[j].indType[0][0] == 'core'):
            jTerm = ops[j]
            jType = jTerm.indices[0].indType
#            if isinstance(iTerm, desOp) and isinstance(jTerm, creOp):
            if isinstance(iTerm, creOp) and isinstance(jTerm, desOp):
#                 print 'kterm=',iTerm,jTerm,iType[0][0],jType[0][0]
#                 exit()
                 if iType[0][0] == 'core' or jType[0][0] == 'core':
##koushik
                    contractionPairs.append((i,j));
#    print "contractionPairs\n", contractionPairs

    # Determine maximum contraction order
    creCount = 0
    maxConOrder = 0
    for i in range(len(ops)-1,-1,-1):
      iTerm = ops[i]
#      if isinstance(iTerm, creOp):
      if isinstance(iTerm, desOp):
        creCount +=1
#      elif isinstance(iTerm, desOp) and creCount > 0:
      elif isinstance(iTerm, creOp) and creCount > 0:
        maxConOrder += 1
        creCount -= 1
    del(creCount,iTerm)

    # Generate all contractions
    contractions = []
    for i in range(maxConOrder+1):
      subCons = makeTuples(i,contractionPairs)
      j = 0
      while j < len(subCons):
        creOpTags = []
        desOpTags = []
        for k in range(i):
          creOpTags.append(subCons[j][k][1])
          desOpTags.append(subCons[j][k][0])
#          creOpTags.append(subCons[j][k][0])
#          desOpTags.append(subCons[j][k][1])
        if allDifferent(creOpTags) and allDifferent(desOpTags):
          j += 1
        else:
          del(subCons[j])
      for j in range(len(subCons)):
        contractions.append(subCons[j])
    del(subCons,creOpTags,desOpTags,contractionPairs)
#    print "contractions:\n", contractions

    # For each contraction, generate the resulting term
    outTerms = []
    for contraction in contractions:
      conSign = 1
      deltaFuncs = []
      conIndeces = []
      subOpString = []
      subOpString.extend(ops)
      for conPair in contraction:
        index1 = ops[conPair[0]].indices[0]
        index2 = ops[conPair[1]].indices[0]
        deltaFuncs.append(kroneckerDelta([index1,index2]))
        subOpString[conPair[0]] = 'contracted'
        subOpString[conPair[1]] = 'contracted'
        for q in subOpString[conPair[0]+1:conPair[1]]:
          if not (q is 'contracted'):
            conSign *= -1
      i = 0
      while i < len(subOpString):
        if subOpString[i] is 'contracted':
          del(subOpString[i])
        else:
          i += 1
      (sortSign,sortedOps) = sortOpsCore(subOpString)
      totalSign = conSign * sortSign
      outTensors = []
      outTensors.extend(nonOps)
      outTensors.extend(deltaFuncs)
      outTensors.extend(sortedOps)
      outTerms.append( term(totalSign * inTerm.numConstant, inTerm.constants, outTensors) )


  # Normal ordering for sfExOps
  elif has_sfExOps:

    # Make separate lists of the spin free excitation operators and other tensors
    raise Exception('This code does not support for now')

  else:
    outTerms = []
    outTerms = [inTerm]
#    return outTerms
#    raise RuntimeError, "Normal ordering function failed to choose what to do."

  print "Terms after normal ordering:"
  for t in outTerms:
    print t

  return outTerms
#
#####################################
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
def sortOpsCore(unsortedOps, returnPermutation = False):
  """
  Sorts a list of creation/destruction operators into normal order and alphebetically.
  Performs no contractions.  Returns the overall sign resulting from the sort and the sorted operator list.
  """
  sortedOps = unsortedOps + []
  i = 0
  sign = 1
#  print 'kperm=',len(unsortedOps),len(sortedOps)
#  for i in sortedOps:
#    print i
#  exit()
#  print 'ki=',i
  if returnPermutation:
    perm = range(len(unsortedOps))
#  print 'kperm=',len(unsortedOps),len(sortedOps)
  while i < len(sortedOps)-1:
####
#    if (sortedOps[i].name == 'cre') and (sortedOps[i].indices[0].indType[0][0]=='core'):
    if (sortedOps[i].name == 'cre') and (sortedOps[i].indices[0].indType[0][0]=='core'):
   #    if sortedOps[i] <= sortedOps[i+1]:
     #     if (sortedOps[i+1].name == sortedOps[i].name):
##########
          j = i
          for k in range(i,(len(sortedOps)-1)):
            if (sortedOps[k+1].name == 'cre') and (sortedOps[k+1].indices[0].indType[0][0]=='core'):
               i += 1
               j = i
            else:
               i = j 
               break
          if ((i+1) > (len(sortedOps)-1)):
              break     
#
#          if (sortedOps[i+1].name == 'cre') and (sortedOps[i+1].indices[0].indType[0][0]=='core'):
##            print 'km0=',sortedOps[i],sortedOps[i+1],i,i+1,len(sortedOps)
#            i +=1
#            if ((i+1) > (len(sortedOps)-1)):
#               break
##########
          temp = sortedOps[i]
          sortedOps[i] = sortedOps[i+1]
          sortedOps[i+1] = temp
#          print 'km1=',sortedOps[i],sortedOps[i+1],i
          if returnPermutation:
            temp = perm[i]
            perm[i] = perm[i+1]
            perm[i+1] = temp
          i = 0
          sign *= -1
     #     if (sortedOps[i+1].name == sortedOps[i].name):
     #       break
    elif (sortedOps[i+1].name == 'des') and (sortedOps[i+1].indices[0].indType[0][0]=='core'):
          temp = sortedOps[i]
          sortedOps[i] = sortedOps[i+1]
          sortedOps[i+1] = temp
#          print 'km2=',sortedOps[i-1],sortedOps[i],i
          if returnPermutation:
            temp = perm[i+1]
            perm[i+1] = perm[i]
            perm[i] = temp
          i = 0
          sign *= -1
          if (sortedOps[i+1].name == sortedOps[i].name):
            break
    else:
#       print "koushik check3",i
       i += 1
  if returnPermutation:
    return (sign,sortedOps,perm)
  return (sign,sortedOps)
#####################################
def print_header():

    print("""\n--------------------------------------------------------------
    SQA_extra: Code geneator for quasi-particle systems.
    author:  Koushik Chatterjee
    date:  August 31, 2018

    Copyright (C) 2018-2019  Koushik Chatterjee (koushikchatterjee7@gmail.com)

    This program is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU General Public License
    for more details.
--------------------------------------------------------------""")
#
#
def contractDeltaFuncs_nondummy(terms):
#
 "Contracts delta function for both non-dummy indices only wrt to orbitals subspaces, otherwise use 'contractDeltaFuncs' function."
#
 for term in terms:
#    i = 0
# while i < len(term.tensors):
    for i in range(len(term.tensors)):
        t = term.tensors[i]
#
        if isinstance(t, kroneckerDelta):
          i0 = t.indices[0]
          i1 = t.indices[1]
          if not (i0.isSummed and i1.isSummed):
              if not (i0.indType[0][0] == i1.indType[0][0]):
                   term.numConstant = 0.0
#
 termChop(terms)
#
 return terms
#
#####################################
def tensorIndex_order(terms):
 print ""
 print "Reorder tensor indices according to (Core < Active < Virtual): =>"
 braSym = True
 ketSym = True
 for term in terms:
     for i in range(len(term.tensors)):
         t = term.tensors[i]
         t0 = t.copy()
         if ((len(t.indices)==2) or (len(t.indices)==4)):
#
            if (len(t.symmetries) == 2):
               if (t.symmetries[0].factor != -1):
                  braSym = False
               if (t.symmetries[1].factor != -1):
                  ketSym = False
#
            ind_list = []
            ind_rank = []
          #  print  term
            for j in range(len(t.indices)):
                ind_list.append(t.indices[j].name)
                if (t.indices[j].indType[0][0][0] == 'c'):
                   rank = 0
                elif(t.indices[j].indType[0][0][0] == 'a'):
                   rank = 1
                else:
                   rank = 3
                ind_rank.append(rank)
            if (len(t.indices)==2):
               if (ind_rank[0] > ind_rank[1]):
                  print  term
                  index_type = t.indices[0].indType
                  t.indices[0].name = ind_list[1]
                  t.indices[0].indType = t.indices[1].indType
                  t.indices[1].name = ind_list[0]
                  t.indices[1].indType = index_type
                  print t0,'    --->', t
            else:
               if (ind_rank[0] > ind_rank[1]):
                  print  term
                  index_type = t.indices[0].indType
                  t.indices[0].name = ind_list[1]
                  t.indices[0].indType = t.indices[1].indType
                  t.indices[1].name = ind_list[0]
                  t.indices[1].indType = index_type
                  if (braSym):
                     term.scale(-1.0)
                  print t0,'    --->', t
               if (ind_rank[2] > ind_rank[3]):
                  index_type = t.indices[2].indType
                  t.indices[2].name = ind_list[3]
                  t.indices[2].indType = t.indices[3].indType
                  t.indices[3].name = ind_list[2]
                  t.indices[3].indType = index_type
                  if (ketSym):
                     term.scale(-1.0)
                  print t0,'    --->', t
               if ((ind_rank[0]+ind_rank[1]) > (ind_rank[2]+ind_rank[3])):
                  print  term
                  t.indices[0].name = ind_list[2]
                  t.indices[0].indType = t.indices[2].indType
                  t.indices[1].name = ind_list[3]
                  t.indices[1].indType = t.indices[3].indType
                  t.indices[2].name = ind_list[0]
                  t.indices[2].indType = t.indices[0].indType
                  t.indices[3].name = ind_list[1]
                  t.indices[3].indType = t.indices[1].indType
                  print t0,'Bra <---> Ket', t

 return terms

#####################################
