#    file:  sqaAddon.py
#  author:  Koushik Chatterjee
#    date:  August 31, 2018
# summary:  Defines the dummy indices according to core, active, external labels. 
#
# (c) 2018-2019 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#

#import sqa_extra.secondQuantizationAlgebra as sqa
from sqaIndex import index
from sqaTensor import tensor, kroneckerDelta, sfExOp, creOp, desOp, creDesTensor
from sqaTerm import term, multiplyTerms, termChop
from sqaOptions import options

# Create an object of the index class
#addon = index()

#####################################
def addon(nterms):
 "A class for dummy indices label upate."
 print "koushik Call Addon"
# coreInd = list('ijklmn')
# actvInd = list('xyzwvz')
# virtInd = list('abcdef')
#
# print coreInd.pop(0)
# print nterms[-1]
# print nterms[-1].tensors[0].indices[2].name


#  Make a dummy list
# ndummy=[]
 bestmap={}
# mymap ={}
##
 for t in nterms:
    ndummy=[]
#    bestmap={}
    mymap ={}
    index_types = ()
#
    coreInd = list('ijklmn')
    actvInd = list('xyzwvz')
    virtInd = list('abcdef')
#
    for t_tensor in t.tensors:
#
#        for t_tensor_index in t_tensor.indices:
#
        for t_tensor_index in range(len(t_tensor.indices)):
#            print "koushik check len=", len(t_tensor.indices), t_tensor.indices[t_tensor_index].name, t_tensor.indices[t_tensor_index].isSummed, t_tensor,t_tensor.indices[t_tensor_index].tup()
#
       #     if t_tensor.indices[t_tensor_index].isSummed == True and t_tensor.indices[t_tensor_index].name not in ndummy:
#
#
            
            #bestmap[t_tensor.indices[t_tensor_index].tup()] = index(t_tensor.indices[t_tensor_index].name, t_tensor.indices[t_tensor_index].indType, t_tensor.indices[t_tensor_index].isSummed)

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
#            print 'ko', t_tensor.indices[t_tensor_index].name,t_tensor.indices[t_tensor_index].indType, t_tensor.indices[t_tensor_index].isSummed, t_tensor.indices[t_tensor_index].tup()
#
#
#                     if t_tensor.indices[t_tensor_index].tup() in bestmap.keys():
#                       bestmap[t_tensor.indices[t_tensor_index].tup()].name = coreInd.pop(0)
#
#            print 'kchk=',t_tensor.indices[t_tensor_index].name
#            print bestmap
#            print 'kchk1=',t_tensor.indices[t_tensor_index].name
#            print bestmap.keys()



#            for key in bestmap.keys():
#               print 'koushik check=',key,t_tensor.indices[t_tensor_index].tup()
#
#            if t_tensor.indices[t_tensor_index].isSummed == True:
#                bestmap[t_tensor.indices[t_tensor_index].tup()].name = coreInd.pop(0)
#                t_tensor.indices[t_tensor_index] = bestmap[t_tensor.indices[t_tensor_index].tup()].copy()
#####################################
##check
    # if t_tensor_index.indType[0] == ('core',):
##check
        #    print bestmap.keys(),t_tensor.indices[t_tensor_index].tup()
#            if bestmap[t_tensor.indices[t_tensor_index].tup()].isSummed == True:
#                for key1 in bestmap.keys():
#                   print key1, bestmap[key1].tup(), bestmap[key1].name, bestmap(1)
#                   return
            #       for key2 in bestmap.keys():
            #          if bestmap[key2].name = bestmap[key1].name:
            #             bestmap[key1].name = coreInd.pop(0) 
         #       bestmap[t_tensor.indices[t_tensor_index].tup()].name = coreInd.pop(0)
         #       t_tensor.indices[t_tensor_index] = bestmap[t_tensor.indices[t_tensor_index].tup()].copy()

#
#
#
#        index_types += t_tensor.indices[t_tensor_index].indType
    print t
#    print index_types
#    print ndummy
#    print "koushik dictionary"
#
#    print 'kbestma=', bestmap.keys()
#####################################
#
#
#            for key in bestmap.keys():
#    #            print 'koushik check=',key,bestmap[key].name,bestmap[key].isSummed
#                if bestmap[key].isSummed == True:
#                   bestmap[key].name = coreInd.pop(0)
#                 #  if key in bestmap.keys():
#                 #    Print "same key" return
#                t_tensor.indices[t_tensor_index] = bestmap[key].copy






#####################################
#              while bestmap.keys():
#                minval = min(bestmap.values())
#                for key in bestmap.keys():
#                   print 'key=',key
#                   if bestmap[key] == minval:
#                     mymap[bestmap[key].tup()] = bestmap[key].copy()
#                     mymap[bestmap[key].tup()].name = coreInd.pop(0)
#           #          print 'kc=',mymap.keys(),mymap[bestmap[key].tup()].name
#                     del bestmap[key]
#                     break
#           #   print 'kc1=',mymap.keys()
#           #   if t_tensor.indices[t_tensor_index].tup() in mymap.keys():
#           #   t_tensor.indices[t_tensor_index] = mymap[t_tensor.indices[t_tensor_index].tup()].copy
                   #        bestmap[t_tensor.indices[t_tensor_index].tup()].name = coreInd.pop(0)
                   #        t_tensor.indices[t_tensor_index] = bestmap[t_tensor.indices[t_tensor_index].tup()].copy()
#####################################
