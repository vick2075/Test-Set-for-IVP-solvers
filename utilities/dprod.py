import numpy as np

## dot product of vectors with unrolled loops

def ddot(a,b,n):
    aa = a.T
    m = np.shape(a)[0]
    if n == m:
        dp = np.dot(aa, b)
    elif n < m and n >= 0:
        s = b[n:]
        s = np.concatenate([s,a[:n,0]])
        dp = np.dot(aa, s)
    else:
        dp = "Error in input"
    return dp




'''
a = np.array([[1, 2, 3, 4],
              [5, 6, 7, 8],
              [9, 10, 11, 12],
              [13, 14, 15, 16]])

b = np.array([17, 18, 19, 20])



dp = np.dot(a, b)
dpt = np.dot(a.T, b)

print("A")
print(a)
print()
print("B")
print(b)
print()
print("Dot Product")
print(dp)
print()
print("Dot Product Column-wise")
print(dpt)
print()

sx = b[3:]
print(sx)
sx = np.concatenate([sx,a[:3,0]])
print(sx)
print()

print("Dot Product Column-wise with unrolled loops")
dptu = np.dot(a.T, sx)
print(dptu)
print()
print('---------------------------------------------------------------------')
print()


print("Algorithm:")
ss = ddot(a, b, 4)
print(ss)
'''
