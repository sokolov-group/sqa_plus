from sqaIndex import index
from sqaTensor import tensor, creOp, desOp, kroneckerDelta
from sqaSymmetry import symmetry
from sqaOptions import options

class creDesTensor(tensor):
  
    freelyCommutes = False 

    def __init__(self, ops, trans_rdm = False):
      
        TypeErrorMessage = "ops must be a normal ordered list of creOp and desOp objects"
        if not type(ops) == type([]):
            raise TypeError, TypeErrorMessage
  
        # Initialize list of cre/des operators
        self.ops = ops

        # Initialize name
        self.trans_rdm = trans_rdm
  
        # Initialize list of cre/des operators
        self.ops = ops

        # Initialize permutations and factors
        (self.permutations, self.factors) = (None, None)
  
        # Count the number of creation/destruction operators
        self.nCre = 0
        self.nDes = 0

        # Build the index list
        self.indices = []
        desFlag = False

        for op in ops:

            # Ensure normal-ordering
            if (not isinstance(op, creOp)) and (not isinstance(op, desOp)):
                raise TypeError, TypeErrorMessage

            if isinstance(op, desOp):
                self.nDes += 1
                desFlag = True

            if isinstance(op, creOp):
                self.nCre += 1
                if desFlag:
                    raise TypeError, TypeErrorMessage

            self.indices.append(op.indices[0].copy())

        # Create count of total indices
        self.nInd = self.nCre + self.nDes

        # Initialize symmetries
        self.symmetries = []
        if len(self.indices) > 1:
            swapValues = range(len(self.indices)-1)
            
            if self.nCre > 0:
                del(swapValues[self.nCre-1])
            
            for i in swapValues:
                if i == 0:
                  temp_tup = (1,)
                else:
                  temp_tup = (0,)
    
                for j in range(1,len(self.indices)):
                    if j == i:
                      temp_tup = temp_tup + (i+1,)
                    elif j == i+1:
                      temp_tup = temp_tup + (i,)
                    else:
                      temp_tup = temp_tup + (j,)
    
                self.symmetries.append(symmetry(temp_tup, -1))

        # Add bra/ket symmetries for ground-state RDMs
        if (len(self.indices) % 2 == 0) and self.trans_rdm == False:
            reversed_range = tuple(range(len(self.indices))[::-1])
            self.symmetries.append(symmetry(reversed_range, 1))

        # Print warning if number of indices is odd and trans_rdm is False
        if (len(self.indices) % 2 != 0) and self.trans_rdm == False:
            print ('trans_rdm flag is set to True, but an ODD number of cre/des operators are present. Switching trans_rdm flag to TRUE !!')
            self.trans_rdm == True

        # Initialize name
        if trans_rdm:
            self.name = 'trdm'
        else:
            self.name = 'rdm'
  

    def __cmp__(self,other):
  
        # comparison to another creDesTensor
        if isinstance(other, creDesTensor):
            retval = cmp(self.name,other.name)
            if retval != 0:
                return retval
            retval = cmp(self.indices,other.indices)
            if retval != 0:
                return retval
            retval = cmp(self.symmetries,other.symmetries)
            return retval
    
        # creDesTensor class is less than the creOp and desOp classes
        elif isinstance(other,creOp):
            return -1
        elif isinstance(other,desOp):
            return -1
    
        # creDesTensor class is greater than other tensor subclasses
        elif isinstance(other,tensor):
            return 1
    
        # Raise an error if other is not a tensor
        else:
            raise TypeError, "A creDesTensor object may only be compared to another tensor"
        return 0
 
 
    def copy(self):

        ops = []

        for i in range(self.nCre):
            ops.append(creOp(self.indices[i]))

        for i in range(self.nCre,len(self.indices)):
            ops.append(desOp(self.indices[i]))

        return creDesTensor(list(ops), self.trans_rdm)
