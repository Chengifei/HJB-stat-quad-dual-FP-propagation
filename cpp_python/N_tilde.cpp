#include "pde_info.h"
#include <span>
using namespace nonlinearity;

template <typename T>
auto get_x_p(const Eigen::PlainObjectBase<T>& args) noexcept {
    return std::make_pair(args.template topRows<M_dim>(),
                          args.template middleRows<L_dim>(M_dim));
}

template <typename T>
Eigen::Matrix<double, 1, T::ColsAtCompileTime> Ntilde::eval(const Eigen::PlainObjectBase<T>& args) {
    auto [x, p] = get_x_p(args);
    return l(x) + (p.array() * f(x).array()).matrix();
}
template <typename T>
Eigen::Matrix<double, mu_dim, T::ColsAtCompileTime> Ntilde::gradient(const Eigen::PlainObjectBase<T>& args) {
    auto [x, p] = get_x_p(args);
    Eigen::Matrix<double, mu_dim, T::ColsAtCompileTime> ret(mu_dim, args.cols());
    ret.template bottomRows<L_dim>() = f(x);
    ret.template topRows<M_dim>() = l_der(x) + (p.array() * f_der(x).array()).matrix();
    return ret;
}
template <typename T, typename U>
auto dot(const T& p, const U& q) noexcept {
    static_assert(T::RowsAtCompileTime == q.size() /* >= 0 */);
    if constexpr (T::RowsAtCompileTime > 1) {
        std::span<typename U::value_type, q.size() - 1> span {q.cbegin() + 1, q.cend()};
        return p[0] * q[0] + dot(p.template middleRows<T::RowsAtCompileTime - 1>(1) - q);
    }
    else
        return p[0] * q[0];
}

Eigen::Matrix<double, mu_dim, mu_dim> Ntilde::dder(const Eigen::Matrix<double, mu_dim, 1>& args) {
    Eigen::Matrix<double, mu_dim, mu_dim> ret;
    auto [x, p] = get_x_p(args);
    ret.bottomRightCorner<L_dim, L_dim>().setZero();
    ret.topLeftCorner<M_dim, M_dim>() = l_dder(x) + dot(p, f_dder(x));
    ret.bottomLeftCorner<L_dim, M_dim>() = f_der(x);
    ret.topRightCorner<M_dim, L_dim>() = ret.bottomLeftCorner<L_dim, M_dim>().transpose();
    return ret;
}

template Eigen::Matrix<double, 1, 1> Ntilde::eval(const Eigen::PlainObjectBase<Eigen::Matrix<double, mu_dim, 1>>&);
template Eigen::Matrix<double, 1, Eigen::Dynamic> Ntilde::eval(const Eigen::PlainObjectBase<Eigen::Matrix<double, mu_dim, Eigen::Dynamic>>&);

template Eigen::Matrix<double, mu_dim, 1> Ntilde::gradient(const Eigen::PlainObjectBase<Eigen::Matrix<double, mu_dim, 1>>&);
template Eigen::Matrix<double, mu_dim, Eigen::Dynamic> Ntilde::gradient(const Eigen::PlainObjectBase<Eigen::Matrix<double, mu_dim, Eigen::Dynamic>>&);
