import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

m = 50
n = 50
a = 15
b = 4
x = np.linspace(0, a, n+1)
y = np.linspace(0, b, m+1)
h1 = a/(n+1)
h2 = b/(m+1)
gamma = 2/h1**2 + 2/h2**2
alpha = (h1/h2)**2
phi = np.zeros((m+1, n+1))

# boundaries
phi[:, 0] = 0
phi[m,:] = np.sin(x)/np.sin(a)
phi[0, :] = 0
phi[:,n] = np.sinh(y)/np.sinh(b)

U = np.copy(phi)

for k in range(10000):
    for i in range(1, m):
        for j in range(1, n):
           U[i, j] = ((U[i-1, j]+U[i+1, j])/h1**2+(U[i, j-1]+U[i, j+1])/h2**2)/gamma

X, Y = np.meshgrid(x, y)
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, Y, U, rstride=1, cstride=1)
plt.show()



k = (n-1)*(m-1)
A = np.zeros((k, k))
N = n
I = np.identity(N - 1)
B = np.zeros((N - 1, N - 1))
np.fill_diagonal(B, -2 * (1 + alpha))
for i in range(N - 2):
    B[i + 1, i] = B[i, i + 1] = alpha
for i in range(0, k, N - 1):
    A[i:i + N - 1, i:i + N - 1] = B

for i in range(0, k - N, N - 1):
    A[i + N - 1:i + N - 2 + N, i:i + N - 1] = A[i:i + N - 1, i + N - 1:i + N - 2 + N] = I

f = np.zeros((m+1, n+1))

F = np.zeros(k)
t = 0

for i in range(1, m):
    for j in range(1, n):
        F[t] = h1 ** 2 * f[i, j]
        if i-1 == 0:
            F[t] -= phi[i - 1, j]
        elif i + 1 == m:
            F[t] -= phi[i + 1, j]
        if j + 1 == n:
            F[t] -= alpha * phi[i, j + 1]
        elif j-1 == 0:
            F[t] -= alpha * phi[i, j - 1]
        t += 1

x = np.linalg.solve(A, F)
phi[1:m, 1:n] = x.reshape(m-1, n-1)

# plt.contourf(X, Y, np.transpose(phi), colorinterpolation, cmap=colourMap)
# ax.plot_surface(X, Y, phi, rstride=1, cstride=1)
# plt.show()