<<<<<<< HEAD
import numpy as np
import time
=======
>>>>>>> 93bf168a30a9d295acb8b2c26984180c2c85dbd0
def matrix_multiply(A, B):
    #dimension check
    if len(A[0]) != len(B):
        raise ValueError("A rows dont match B columns")
    #stores results
    multiplied = []
    for i in range(len(A)):
        row = []
        for j in range(len(B[0])):
            product = 0
            #multiplies and sums products
            for k in range(len(B)):
                product += A[i][k] * B[k][j]
            row.append(product)
        #stores multiplication results
        multiplied.append(row)
    return multiplied



A = [[1,2], [3,4],[6,7]]
B = [[5,6], [7,8]]
<<<<<<< HEAD
start_time = time.perf_counter()
func_product = matrix_multiply(A,B)
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(func_product)
print(f"Function elpased time {elapsed_time}")



#2.2 Product using numpy
np_product = np.dot(A,B)

start_time = time.perf_counter()
np_product = np.dot(A,B)
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(np_product)
print(f"Numpy elpased time {elapsed_time}")

=======
print(matrix_multiply(A,B))


>>>>>>> 93bf168a30a9d295acb8b2c26984180c2c85dbd0
