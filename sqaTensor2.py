from sqaIndex import index
from sqaTensor import tensor, creOp, desOp, kroneckerDelta
from sqaSymmetry import symmetry

class creDesTensor(tensor):
  
    freelyCommutes = False
    name = "creDesTensor"
  
    def __init__(self, name, ops, transRDM = False):
      
        TypeErrorMessage = "ops must be a normal ordered list of creOp and desOp objects"
        if not type(ops) == type([]):
            raise TypeError, TypeErrorMessage
  
        # Initialize name
        self.name = name
  
        # Initialize permutations and factors
        (self.permutations, self.factors) = (None, None)
  
        # Build the index list and count the number of creation operators
        self.nCre = 0
        self.indices = []
        desFlag = False

        for op in ops:

            if (not isinstance(op, creOp)) and (not isinstance(op, desOp)):
                raise TypeError, TypeErrorMessage

            if isinstance(op, desOp):
                desFlag = True

            if isinstance(op, creOp):
                self.nCre += 1
                if desFlag:
                    raise TypeError, TypeErrorMessage

            self.indices.append(op.indices[0].copy())
  
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

        if len(self.indices) == 4:
            self.symmetries.append(symmetry((2,3,0,1), 1))

        if len(self.indices) == 6:
                self.symmetries.append(symmetry((3,4,5,0,1,2), 1))
   
 
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

        return creDesTensor(self.name, ops)
