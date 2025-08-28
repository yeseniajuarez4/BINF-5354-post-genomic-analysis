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
print(matrix_multiply(A,B))


