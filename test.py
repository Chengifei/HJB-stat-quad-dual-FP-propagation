import HJB_solve
import numpy
import itertools
import concurrent.futures
import math

def x(a, b):
    ax = numpy.linspace(-(l - 1) / 2., (l - 1) / 2., l)
    gauss = math.sqrt(0.5 / math.pi) * numpy.exp(-0.5 * numpy.square(ax))
    kernel = b * numpy.outer(gauss, gauss)
    return numpy.require((kernel + a).flatten(), 'd', 'A')

l = 5
N = (l * l) // 2
A = numpy.identity(l * l) * -2
A += numpy.diagflat(numpy.full((l * l - 1, ), 0.5), 1)
A += numpy.diagflat(numpy.full((l * l - 1, ), 0.5), -1)
A += numpy.diagflat(numpy.full((l * l - l, ), 0.5), l)
A += numpy.diagflat(numpy.full((l * l - l, ), 0.5), -l)
A = numpy.require(A, 'd', 'FA')
G_x = numpy.zeros((l * l, l * l), 'd')
for r in range(l):
    G_x[r, r + 1] = 1
    G_x[r, r] = -1
    for c in range(1, l - 1):
        G_x[r + l * c, r + l * (c + 1)] = 0.5
        G_x[r + l * c, r + l * (c - 1)] = -0.5
    G_x[r + l * (l - 1), r + l * (l - 1)] = 1
    G_x[r + l * (l - 1), r + l * (l - 2)] = -1
G_y = numpy.zeros((l * l, l * l), 'd')
for c in range(l):
    G_y[l * c, l * c + 1] = 1
    G_y[l * c, l * c] = -1
    for r in range(1, l - 1):
        G_y[r + l * c, r + l * c + 1] = 0.5
        G_y[r + l * c, r + l * c - 1] = -0.5
    G_y[l - 1 + l * c, l - 1 + l * c] = 1
    G_y[l - 1 + l * c, l - 2 + l * c] = -1
C = -0.2 * (G_x.T @ G_x + G_y.T @ G_y)
C = numpy.require(C, 'd', 'FA')
M_0 = x(0, 1).reshape((1, 25))
M_0 = numpy.require(M_0, 'd', 'FA')
L_0 = 0.2 * numpy.ones((25, 1))
C_hat = numpy.asfortranarray([[4, 0], [0, -3]])
C_hat = numpy.require(C_hat, 'd', 'FA')
sigma = numpy.zeros((25, 1))
sigma[0:l] = 1
sigma[-l:] += 1
sigma[0::l] += 1
sigma[(l - 1)::l] += 1
sigma = sigma
Gamma = sigma @ sigma.T
Gamma += C_hat[1, 1] * L_0 @ L_0.T
Gamma = numpy.require(Gamma, 'd', 'FA')
t = numpy.linspace(0, 1, 6)

solver = HJB_solve.OptCtrlProb(A, C, M_0, L_0, C_hat, Gamma)
x1 = numpy.linspace(-10, 10, 51)
x2 = numpy.linspace(-10, 10, 51)
res = numpy.empty((x1.size, x2.size), dtype='d', order='F')
err = numpy.empty((x1.size, x2.size), dtype='d', order='F')
grad = numpy.empty((A.shape[0], x1.size, x2.size), dtype='d', order='F')
pde_err = numpy.empty((x1.size, x2.size), dtype='d', order='F')
N = numpy.asarray([120, 120, 100, 100, 100], 'uint64')

with concurrent.futures.ThreadPoolExecutor(max_workers=15) as executor:
    tasks = {
        executor.submit(solver.solve, x(a, b), t, N): (i, j)
        for (i, a), (j, b) in itertools.product(enumerate(x1), enumerate(x2))
    }
    for future in concurrent.futures.as_completed(tasks):
        i, j = tasks[future]
        try:
            W, e, Wx = future.result()
            print(e)
        except Exception as exc:
            print(f'{i}, {j} generated an exception: {exc}')
        else:
            res[i, j] = W
            err[i, j] = e
            grad[:, i, j] = Wx

with open('res', 'wb') as f:
    f.write(res.tobytes('F'))

with open('err', 'wb') as f:
    f.write(err.tobytes('F'))

with open('grad', 'wb') as f:
    f.write(grad.tobytes('F'))

res2 = numpy.empty((x1.size, x2.size), dtype='d', order='F')

with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
    t2 = t.copy()
    t2[0] = 0.005
    tasks = {
        executor.submit(solver.solve, x(a, b), t2, N): (i, j)
        for (i, a), (j, b) in itertools.product(enumerate(x1), enumerate(x2))
    }
    for future in concurrent.futures.as_completed(tasks):
        i, j = tasks[future]
        try:
            W, e, Wx = future.result()
            print(e)
        except Exception as exc:
            print(f'{i}, {j} generated an exception: {exc}')
        else:
            res2[i, j] = W

Wt = (res2 - res) / t2[0]

for (i, a), (j, b) in itertools.product(enumerate(x1), enumerate(x2)):
    y = x(a, b)
    g = grad[:, i, j]
    terms = [Wt[i, j], 0.5 * y.dot(C @ y), g.dot(A @ y), -0.5 * (sigma.T @ g).dot(sigma.T @ g), HJB_solve.N_tilde(M_0 @ y, L_0.T @ g)[0]]
    print(terms)
    pde_err[i, j] = sum(terms) / sum(abs(i) for i in terms)

with open('pde_err', 'wb') as f:
    f.write(pde_err.tobytes('F'))

from matplotlib import pyplot
fig, (ax1, ax2, ax3) = pyplot.subplots(1, 3, subplot_kw={'projection': '3d'})
X1, X2 = numpy.meshgrid(x1, x2)
ax1.plot_surface(X1, X2, res)
ax1.set_xlabel('x_1')
ax1.set_ylabel('x_2')
cs = ax2.plot_surface(X1, X2, err, cmap='Reds')
fig.colorbar(cs, ax=ax2, shrink=0.5)
ax3.plot_surface(X1, X2, numpy.linalg.norm(grad, axis=0))
fig.savefig("output.png")
pyplot.show()
