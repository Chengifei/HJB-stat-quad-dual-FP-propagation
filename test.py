import HJB_solve
import numpy
import itertools
import concurrent.futures

A = numpy.asfortranarray([[-4, 4, 0], [4, -4, 5], [0, -2, 0]], 'd')
C = numpy.asfortranarray([[0.5, -0.4, 0], [-0.4, 0.5, 0], [0, 0, 0.5]], 'd')
M_0 = numpy.asfortranarray([[1, 0, 0]], 'd')
L_0 = numpy.asfortranarray([[1], [0], [0]], 'd')
C_hat = numpy.asfortranarray([[5.5, 0], [0, -1]], 'd')
Gamma = numpy.asfortranarray([[0, 0, 0], [0, 0, 0], [0, 0, 1]], 'd')
t = numpy.linspace(0, 2, 201)

solver = HJB_solve.OptCtrlProb(A, C, M_0, L_0, C_hat, Gamma, t)
x1 = numpy.linspace(-5, 5, 51)
x2 = numpy.linspace(-5, 5, 51)
res = numpy.empty((51, 51), dtype='d', order='F')
err = numpy.empty((51, 51), dtype='d', order='F')
grad = numpy.empty((3, 51, 51), dtype='d', order='F')
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
    tasks = {
        executor.submit(solver.solve, numpy.asarray((a, b, -1), 'd'), 45): (i, j)
        for (i, a), (j, b) in itertools.product(enumerate(x1), enumerate(x2))
    }
    for future in concurrent.futures.as_completed(tasks):
        i, j = tasks[future]
        try:
            W, e, Wx = future.result()
        except Exception as exc:
            print(f'{i}, {j} generated an exception: {exc}')
        else:
            res[i, j] = W[-1]
            err[i, j] = e[-1]
            grad[:, i, j] = Wx[-1, :]

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
ax2.plot_surface(X1, X2, err, cmap='Reds')
ax3.plot_surface(X1, X2, numpy.linalg.norm(grad, axis=0))
fig.savefig("output.png")
pyplot.show()
