import numpy as np

# comparison is elementwise within np tolerances

class matrixgroupelement(object):
    """
    Matrices as group elements. 
    """
    def __init__(self,M):
        if np.linalg.det(M)==0:
            raise ValueError('input has determinant 0')
        else:
            self.M=M

    def __mul__(self,other):
        return matrixgroupelement(self.M.dot(other.M))

    def __pow__(self,n):
        return matrixgroupelement(np.linalg.matrix_power(self.M,n))

    def __eq__(self,other):
        return np.array_equal(self.M,other.M)

    def __neq__(self,other):
        return not self==other

    def __str__(self):
        return '[['+str(self.M[0])+'],['+str(self.M[1])+']]'

    def __getitem__(self,sliced):
        return self.M[sliced]
    
    def __repr__(self):
        return self.__str__()






if __name__ == "__main__":
    import doctest
    doctest.testmod()
