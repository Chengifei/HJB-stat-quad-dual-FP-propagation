// Copyright (C) 2024 Y. Zheng
// SPDX-License-Identifier: BSD-3-Clause
#include "solver.h"
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
solver::solver(T1&& A, T2&& C, T3&& M_0, T4&& L_0, T5&& C_hat, T6&& Gamma, T7&& t) {
    info.A = A;
    info.C = C;
    info.M_0 = M_0;
    info.L_0 = L_0;
    info.C_hat = C_hat;
    info.C_inv = info.C_hat.fullPivLu();
    new (&info._) Ntilde(-info.C_hat);
    info.Gamma = Gamma;
    info.t = t;
    info.A_bar.topLeftCorner<dim, dim>(info.A.rows(), info.A.cols()) = info.A;
    info.A_bar.bottomLeftCorner<dim, dim>(info.C.rows(), info.C.cols()) = info.C_hat(0, 0) * (info.M_0).transpose() * info.M_0 - info.C;
    info.A_bar.bottomRightCorner<dim, dim>(info.A.cols(), info.A.rows()) = -info.A.transpose();
    info.A_bar.topRightCorner<dim, dim>(info.Gamma.rows(), info.Gamma.cols()) = -info.Gamma;
    info.STM(cache);
    mult.setZero();
    mult.topLeftCorner(info.M_0.rows(), dim) = info.M_0;
    mult.bottomRightCorner(info.L_0.cols(), dim) = info.L_0.transpose();
}
#ifdef MEX
#include "mex.hpp"
#include "mexAdapter.hpp"
class MexFunction : public matlab::mex::Function, solver {
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    matlab::data::ArrayFactory factory;
public:
    MexFunction() : matlabPtr(getEngine()) {}
    static std::vector<double>
        matlab_1_to_eigen(const matlab::data::TypedArray<double>& arr) {
        std::vector<double> ret(arr.getNumberOfElements());
        std::copy_n(arr.cbegin(), ret.size(), ret.data());
        return ret;
    }
    static Eigen::MatrixXd
        matlab_2_to_eigen(const matlab::data::TypedArray<double>& arr) {
        matlab::data::ArrayDimensions d = arr.getDimensions();
        assert(d.size() == 2 && arr.getType() == matlab::data::ArrayType::DOUBLE);
        Eigen::MatrixXd ret(d[0], d[1]);
        std::copy_n(arr.cbegin(), ret.size(), ret.data());
        return ret;
    }
    template <typename T = Eigen::MatrixXd>
    static std::vector<T> matlab_3_to_eigen(const matlab::data::TypedArray<double>& arr) {
        matlab::data::ArrayDimensions d = arr.getDimensions();
        std::vector<T> ret;
        ret.reserve(d[2]);
        auto it = arr.cbegin();
        for (std::size_t i = 0; i != d[2]; ++i) {
            Eigen::MatrixXd tmp(d[0], d[1]);
            std::copy_n(it, tmp.size(), tmp.data());
            std::advance(it, tmp.size());
            ret.push_back(std::move(tmp));
        }
        return ret;
    }
    template <typename T = Eigen::MatrixXd>
    static std::vector<std::vector<T>>
        matlab_4_to_eigen(const matlab::data::TypedArray<double>& arr) {
        matlab::data::ArrayDimensions d = arr.getDimensions();
        std::vector<std::vector<T>> ret;
        ret.reserve(d[3]);
        auto it = arr.cbegin();
        for (std::size_t i = 0; i != d[3]; ++i) {
            std::vector<T> tmp0;
            tmp0.reserve(d[2]);
            for (std::size_t j = 0; j != d[2]; ++j) {
                Eigen::MatrixXd tmp(d[0], d[1]);
                std::copy_n(it, tmp.size(), tmp.data());
                std::advance(it, tmp.size());
                tmp0.push_back(std::move(tmp));
            }
            ret.push_back(std::move(tmp0));
        }
        return ret;
    }
    template <typename T>
    auto eigen_to_matlab(Eigen::PlainObjectBase<T>& ret) {
        return factory.createArray({ static_cast<std::size_t>(ret.rows()), static_cast<std::size_t>(ret.cols()) }, ret.data(), ret.data() + ret.size());
    }
    void construct(const matlab::data::Array& _info, const matlab::data::Array& C_hat,
            const matlab::data::Array& D_bar, const matlab::data::Array& A_bar) {
        //verify(t.getDimensions().size() <= 2 && t.getType() == matlab::data::ArrayType::DOUBLE);
        verify(D_bar.getDimensions().size() == 3 && D_bar.getType() == matlab::data::ArrayType::DOUBLE);
        verify(A_bar.getDimensions().size() == 3 && A_bar.getType() == matlab::data::ArrayType::DOUBLE);
        static_cast<solver*>(this)->~solver();
        new (static_cast<solver*>(this)) solver(
            matlab_2_to_eigen(matlabPtr->getProperty(_info, "A")),
            matlab_2_to_eigen(matlabPtr->getProperty(_info, "C")),
            matlab_2_to_eigen(matlabPtr->getProperty(_info, "M0")),
            matlab_2_to_eigen(matlabPtr->getProperty(_info, "L0")),
            matlab_2_to_eigen(C_hat),
            matlab_2_to_eigen(matlabPtr->getProperty(_info, "Gamma")),
            matlab_1_to_eigen(matlabPtr->getProperty(_info, "t")));
    }
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        if (inputs.empty())
            matlabPtr->feval(u"error",
                0, std::vector<matlab::data::Array>({ factory.createScalar("Too few arguments: No command") }));
        if (inputs[0].getType() != matlab::data::ArrayType::MATLAB_STRING)
            matlabPtr->feval(u"error",
                0, std::vector<matlab::data::Array>({ factory.createScalar("The first argument <command> must be a char array") }));
        matlab::data::TypedArray<matlab::data::MATLABString> _cmd(inputs[0]);
        std::string cmd = _cmd[0];
        if (cmd == "new") [[unlikely]] {
            verify(inputs.size() == 5);
            construct(inputs[1], inputs[2], inputs[3], inputs[4]);
        }
        else if (cmd == "z_Pz+q") {
            verify(inputs.size() == 4);
            auto ret = qint_LTI(matlab_2_to_eigen(inputs[1]), matlab_2_to_eigen(inputs[2]), info.t.cbegin() + static_cast<Eigen::Index>(inputs[3][0]), info.t.cbegin(), info.t.cend());
            outputs[0] = eigen_to_matlab(ret);
        }
        else if (cmd == "N_tilde") {
            verify(inputs.size() == 2);
            if (outputs.empty())
                return;
            auto ret = info._.eval(matlab_2_to_eigen(inputs[1]));
            outputs[0] = eigen_to_matlab(ret);
        }
        else if (cmd == "eta_inv") {
            verify(inputs.size() == 2);
            if (outputs.empty())
                return;
            auto ret = info.eta_inv(matlab_2_to_eigen(inputs[1]));
            outputs[0] = factory.createArray({ static_cast<std::size_t>(ret.rows()), static_cast<std::size_t>(ret.cols()) }, ret.data(), ret.data() + ret.size());

        }
        else if (cmd == "W") {
            verify(inputs.size() == 3);
            if (outputs.empty())
                return;
            try {
                auto [z_Pzq, free_resp] = solve(matlab_2_to_eigen(inputs[1]), static_cast<Eigen::Index>(inputs[2][0]));
                auto ret = info.eta_inv(mult * z_Pzq);
                outputs[0] = factory.createArray({ static_cast<std::size_t>(ret.rows()), static_cast<std::size_t>(ret.cols()) }, ret.data(), ret.data() + ret.size());
                if (outputs.size() > 1)
                    outputs[1] = eigen_to_matlab(free_resp);
            }
            catch (iterative_solver_exc& err) {
                char buf[64];
                std::snprintf(buf, sizeof buf, "%s fp_exc: i = %llu, e = %e", err.str, err.i, err.e);
                matlabPtr->feval(u"error",
                    0, std::vector<matlab::data::Array>({ factory.createScalar(buf) }));
            }
        }
        else
            matlabPtr->feval(u"error",
                0, std::vector<matlab::data::Array>({ factory.createScalar("Unknown command") }));
    }
};
#else
#include "Python.h"
#include "numpy/arrayobject.h"

struct _solver : PyObject, solver {};

static PyModuleDef Module = {
    PyModuleDef_HEAD_INIT,
    "HJB_solve",
    "",
    -1,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr
};

void dealloc(PyObject* self) noexcept {
    static_cast<_solver*>(self)->~_solver();
}

template <typename T>
T numpy_to_eigen(PyArrayObject* arr) {
    npy_intp* dim = PyArray_SHAPE(arr);
    return Eigen::Map<T>(static_cast<T::Scalar*>(PyArray_DATA(arr)), dim[0], dim[1]);
}

struct Python_API_Exception {};

template <class T, class U>
T PyExc(T t, U fail_ret) {
    if (t == fail_ret)
        throw Python_API_Exception();
    return t;
}

template <class T, class U>
T PyOnly(T t, U good_ret) {
    if (t != good_ret)
        throw Python_API_Exception();
    return t;
}

int init(PyObject* self, PyObject* args, PyObject* kwds) noexcept {
    static const char* keys[] = { "A", "C", "M_0", "L_0", "C_hat", "Gamma", "t", nullptr };
    PyArrayObject* arr[7];
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOOOO", const_cast<char**>(keys), arr, arr + 1, arr + 2, arr + 3, arr + 4, arr + 5, arr + 6))
        return -1;
    for (auto i : arr) {
        if (!PyArray_Check(i)) {
            PyErr_BadArgument();
            return -1;
        }
        if (!PyArray_IS_F_CONTIGUOUS(i)) {
            PyErr_BadArgument();
            return -1;
        }
        if (PyArray_TYPE(i) != NPY_DOUBLE) {
            PyErr_BadArgument();
            return -1;
        }
    }
    for (int i = 0; i != 6; ++i)
        if (PyArray_NDIM(arr[i]) != 2) {
            PyErr_BadArgument();
            return -1;
        }
    std::vector<double> t(static_cast<double*>(PyArray_DATA(arr[6])), static_cast<double*>(PyArray_DATA(arr[6])) + PyArray_SIZE(arr[6]));
    using namespace Eigen;
    new (static_cast<solver*>(static_cast<_solver*>(self))) solver(
        numpy_to_eigen<Matrix<double, dim, dim, DontAlign>>(arr[0]),
        numpy_to_eigen<Matrix<double, dim, dim, DontAlign>>(arr[1]),
        numpy_to_eigen<Matrix<double, Dynamic, dim, DontAlign>>(arr[2]),
        numpy_to_eigen<Matrix<double, dim, Dynamic, DontAlign>>(arr[3]),
        numpy_to_eigen<Matrix<double, mu_dim, mu_dim, DontAlign>>(arr[4]),
        numpy_to_eigen<Matrix<double, dim, dim, DontAlign>>(arr[5]),
        std::move(t)
    );
    return 0;
}

PyObject* solve(PyObject* self, PyObject* const args[], Py_ssize_t argc) {
    Eigen::Matrix<double, dim, 1> x;
    if (argc < 2)
        return nullptr;
    PyArrayObject* const _x = reinterpret_cast<PyArrayObject*>(args[0]);
    Eigen::Index T = PyLong_AsLong(args[1]);
    if (!PyArray_Check(_x)) {
        PyErr_BadArgument();
        return nullptr;
    }
    if (PyArray_TYPE(_x) != NPY_DOUBLE) {
        PyErr_BadArgument();
        return nullptr;
    }
    std::copy(static_cast<double*>(PyArray_DATA(_x)), static_cast<double*>(PyArray_DATA(_x)) + dim, x.data());
    PyThreadState* _save;
    Py_UNBLOCK_THREADS
    auto [W, err, Wx] = static_cast<_solver*>(self)->W(x, T);
    Py_BLOCK_THREADS
    npy_intp dim = W.size();
    assert(err.size() == W.size());
    PyArrayObject* py1 = reinterpret_cast<PyArrayObject*>(PyExc(PyArray_SimpleNew(1, &dim, NPY_DOUBLE), nullptr));
    std::copy_n(W.data(), W.size(), static_cast<double*>(PyArray_DATA(py1)));
    PyArrayObject* py2 = reinterpret_cast<PyArrayObject*>(PyExc(PyArray_SimpleNew(1, &dim, NPY_DOUBLE), nullptr));
    std::copy_n(err.data(), err.size(), static_cast<double*>(PyArray_DATA(py2)));
    npy_intp dims[2] = { Wx.cols(), ::dim};
    PyArrayObject* py3 = reinterpret_cast<PyArrayObject*>(PyExc(PyArray_SimpleNew(2, dims, NPY_DOUBLE), nullptr));
    std::copy_n(Wx.data(), Wx.size(), static_cast<double*>(PyArray_DATA(py3)));
    auto py_ret = PyTuple_Pack(3, py1, py2, py3);
    Py_DECREF(py1);
    Py_DECREF(py2);
    Py_DECREF(py3);
    return py_ret;
}

static PyMethodDef methods[] = {
    {"solve", reinterpret_cast<PyCFunction>(solve), METH_FASTCALL, "solve(x), solve for the mu process passing through (t[0], x)."},
    {nullptr}
};

PyTypeObject optimal_control_problem_type {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "HJB_solve.OptCtrlProb",
    .tp_basicsize = sizeof(_solver),
    .tp_itemsize = 0,
    .tp_dealloc = dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_methods = methods,
    .tp_init = init,
    .tp_new = PyType_GenericNew
};

PyMODINIT_FUNC
PyInit_HJB_solve() {
    PyObject* m = PyExc(PyModule_Create(&Module), nullptr);
    import_array();
    PyOnly(PyType_Ready(&optimal_control_problem_type), 0);
    if (PyModule_AddObjectRef(m, "OptCtrlProb", reinterpret_cast<PyObject*>(&optimal_control_problem_type))) {
        Py_XDECREF(&optimal_control_problem_type);
        return nullptr;
    }
    return m;
}

#endif
