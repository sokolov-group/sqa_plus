from sqaIndex import index
from sqaTensor import tensor, creOp, desOp, kroneckerDelta
from sqaSymmetry import symmetry
from sqaOptions import options

class creDesTensor(tensor):
  
    freelyCommutes = False
    name = 'rdm'
 
    def __init__(self, name, ops, transRDM = False):
      
        TypeErrorMessage = "ops must be a normal ordered list of creOp and desOp objects"
        if not type(ops) == type([]):
            raise TypeError, TypeErrorMessage
  
        # Initialize name
        self.name = name
  
        # Initialize list of cre/des operators
        self.ops = ops

        # Initialize name
        self.transRDM = transRDM
  
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
        if (len(self.indices) % 2 == 0) and self.transRDM == False:
            reversed_range = tuple(range(len(self.indices))[::-1])
            self.symmetries.append(symmetry(reversed_range, 1))

        # Add bra/ket symmetries for ground-state RDMs
        if (len(self.indices) % 2 != 0) and self.transRDM == False:
            print ('transRDM flag is set to True, but an ODD number of cre/des operators are present. Switching transRDM flag to TRUE !!')
            self.transRDM == True


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

        print (ops)
        return creDesTensor(self.name, list(ops))
