import HJB_solve
import numpy
import itertools
import concurrent.futures

def x(a, b):
    ax = numpy.linspace(-(l - 1) / 2., (l - 1) / 2., l)
    gauss = numpy.exp(-0.5 * numpy.square(ax))
    kernel = (b - a) * numpy.outer(gauss, gauss)
    return numpy.require((kernel + a).flatten(), 'd', 'A')

l = 5
N = (l * l) // 2
A = numpy.identity(l * l) * -2
A += numpy.diagflat(numpy.full((l * l - 1, ), 0.5), 1)
A += numpy.diagflat(numpy.full((l * l - 1, ), 0.5), -1)
A += numpy.diagflat(numpy.full((l * l - l, ), 0.5), l)
A += numpy.diagflat(numpy.full((l * l - l, ), 0.5), -l)
A *= 0.2
A = numpy.require(A, 'd', 'FA')
C = 0.11 * numpy.identity(l * l)
for j in range(l * l):
    if N != j:
        C[N, j] = -0.1
        C[j, N] = -0.1
C[N, N] = 0.1 * l * l + 0.1
C = numpy.require(C, 'd', 'FA')
M_0 = x(0, 1).reshape((1, l * l))
M_0 = numpy.require(M_0, 'd', 'FA')
L_0 = numpy.zeros((l * l, 1))
L_0[N] = 1
C_hat = numpy.asfortranarray([[4, 0], [0, -2]])
C_hat = numpy.require(C_hat, 'd', 'FA')
sigma = numpy.zeros((l * l, 4))
weights = numpy.linspace(0, 1, l + 2)[1:-1]
sigma[:l, 0] = 1 - weights
sigma[:l, 1] = weights
sigma[(l - 1)::l, 1] += 1 - weights
sigma[(l - 1)::l, 2] = weights
sigma[(-l):, 2] += weights
sigma[(-l):, 3] = 1 - weights
sigma[:((l - 1) * l + 1):l, 0] += 1 - weights
sigma[:((l - 1) * l + 1):l, 3] += weights
Gamma = sigma @ sigma.T
Gamma += C_hat[1, 1] * L_0 @ L_0.T
Gamma = numpy.require(C, 'd', 'FA')
t = numpy.linspace(0, 1, 5)

solver = HJB_solve.OptCtrlProb(A, C, M_0, L_0, C_hat, Gamma)
x1 = numpy.linspace(-2, 2, 21)
x2 = numpy.linspace(-2, 2, 21)
res = numpy.empty((x1.size, x2.size), dtype='d', order='F')
err = numpy.empty((x1.size, x2.size), dtype='d', order='F')
grad = numpy.empty((A.shape[0], x1.size, x2.size), dtype='d', order='F')
N = numpy.asarray([120, 100, 100, 50], 'uint64')

with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
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
