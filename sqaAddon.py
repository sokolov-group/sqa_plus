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

import sqa_extra.secondQuantizationAlgebra as sqa
from sqaIndex import index
from sqaTensor import tensor, kroneckerDelta, sfExOp, creOp, desOp, creDesTensor
from sqaTerm import term, multiplyTerms, termChop
from sqaOptions import options

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
