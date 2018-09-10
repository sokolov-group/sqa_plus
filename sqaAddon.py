#    file:  sqaAddon.py
#  author:  Koushik Chatterjee
#    date:  August 31, 2018
# summary:  addon: Defines the dummy indices according to core, active, external labels. 
#          addon1: Evaluate expectation value and kill the zero term (make the coefficient 0). if some one want to delete, may call termChop.
#
# (c) 2018-2019 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#

from sqaIndex import index
from sqaTensor import tensor, kroneckerDelta, sfExOp, creOp, desOp, creDesTensor
from sqaTerm import term, multiplyTerms, termChop, sortOps
from sqaOptions import options
from sqaMisc import makeTuples, allDifferent, makePermutations

from sqaNormalOrder import normalOrder

# Create an object of the index class
#addon = index()

#####################################
def addon(nterms):
 "A function for dummy indices label upate."
#
 print "Dummy indices label update:=>"
##
 for t in nterms:
    mymap ={}
    index_types = ()
#
    coreInd = list('ijklmn')
    actvInd = list('xyzwvz')
    virtInd = list('abcdef')
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
 print ""
 print "A function for  expectation value: kill the zero terms."
 print ""
 addon1(nterms)
#
#####################################
#
def addon1(nterms):
 "A function for  expectation value: kill the zero terms."
#
 print "Expectation value ( Kill the zero terms):=>"
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
def normalOrderCore(inTerm):
  "Returns a list of terms resulting from normal ordering the operators in inTerm."

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
    print "contractions:\n", contractions

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
      (sortSign,sortedOps) = sortOps1(subOpString)
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
    return
    raise RuntimeError, "Normal ordering function failed to choose what to do."

  print "Terms after normal ordering:"
  for t in outTerms:
    print t

  return outTerms


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
def sortOps1(unsortedOps, returnPermutation = False):
  """
  Sorts a list of creation/destruction operators into normal order and alphebetically.
  Performs no contractions.  Returns the overall sign resulting from the sort and the sorted operator list.
  """
  sortedOps = unsortedOps + []
  i = 0
  sign = 1
#  print 'kperm=',len(unsortedOps),len(sortedOps)
  if returnPermutation:
    perm = range(len(unsortedOps))
#    print 'kperm=',perm
  while i < len(sortedOps)-1:
    if sortedOps[i] >= sortedOps[i+1]:
#       print 'km=',sortedOps[i],sortedOps[i+1]
       i += 1
    else:
#      print 'km1=',sortedOps[i],sortedOps[i+1]
      temp = sortedOps[i+1]
      sortedOps[i+1] = sortedOps[i]
      sortedOps[i] = temp
      if returnPermutation:
        temp = perm[i+1]
        perm[i+1] = perm[i]
        perm[i] = temp
      i = 0
      sign *= -1
  if returnPermutation:
    return (sign,sortedOps,perm)
  return (sign,sortedOps)

