#include "solver.h"
#include "transform.h"

std::vector<double> linspace(double begin, double end, std::size_t size) {
    if (size == 0)
        return {};
    std::vector<double> ret(size);
    const double step = (end - begin) / size;
    for (std::size_t i = 0; i != size; ++i)
        ret[i] = begin + step * i;
    ret[size - 1] = end;
    return ret;
}

std::vector<double> logspace(double begin, double end, std::size_t size) {
    auto lin = linspace(begin, end, size);
    std::transform(lin.begin(), lin.end(), lin.begin(), [](double e) { return std::exp2(e); });
    return lin;
}

void PDE_info::STM(cache_t& cache) const {
	cache.LTI_STM.resize(t.size());
		using namespace boost::numeric::odeint;
		auto stepper = make_dense_output<runge_kutta_dopri5<Eigen::Matrix<double, 2 * dim, 2 * dim>>>(1e-10, 1e-12);
		std::size_t counter = 0;
		Eigen::Matrix<double, 2 * dim, 2 * dim> identity = Eigen::Matrix<double, 2 * dim, 2 * dim>::Identity();
		integrate_times(stepper, no_control{ *this }, identity, t.cbegin(), t.cend(), 1e-5,
			[&](const auto& state, double t) {
				cache.LTI_STM[counter++] = state;
			});
}

struct exception {
    std::string t;
};

void verify(bool b, std::string t = {}) {
    if (!b)
        throw exception{ std::move(t) };
}

typedef Eigen::Matrix<double, mu_dim, Eigen::Dynamic> mu_matrix;
std::tuple<double, double, Eigen::Matrix<double, dim, 1>>
solver::W(const Eigen::Matrix<double, dim, 1>& x, std::size_t T) const {
    const Eigen::Index _dim = info.A.rows();
    std::vector<double> errs = logspace(-3, -23, 5);
    std::vector<std::size_t> N(info.t.size() - 1, 50);
    auto [w, grad, err] = solve(x, info.t, N, errs);
    return std::make_tuple(w, err.norm(), grad);
}
template <std::random_access_iterator Iter>
Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>
solver::qint_LTI(const Eigen::Matrix<double, dim, 1>& x, const mu_matrix& mu_array, Iter T, Iter begin) const {
    return qint_LTI(x, mu_array, T, begin, T);
}
template <std::random_access_iterator Iter>
Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>
solver::qint_LTI(const Eigen::Matrix<double, dim, 1>& x, const mu_matrix& mu_array, Iter T, const Iter begin, const Iter end) const {
    using namespace boost::numeric::odeint;
    assert(T >= begin && T <= end);
    Eigen::Matrix<double, 2 * dim, Eigen::Dynamic> ret(2 * info.A.rows(), end - begin);
    typedef Eigen::Matrix<double, 2 * dim, 1> state_t;
    if (T > begin) [[likely]] {
        auto stepper = make_controlled<runge_kutta_dopri5<state_t>>(1e-10, 1e-10);
        state_t zeros = Eigen::Matrix<double, 2 * dim, 1>::Zero(2 * info.A.rows());
        Eigen::Index counter = T - begin;
        integrate_times(stepper, backward{ info, mu_array },
            zeros,
            std::make_reverse_iterator(T), std::make_reverse_iterator(begin), -1e-5,
            [&](const state_t& state, double t) {
                ret.col(--counter) = state;
            });
        }
    if (T < end) {
        auto stepper = make_dense_output<runge_kutta_dopri5<state_t>>(1e-10, 1e-10);
        state_t zeros = Eigen::Matrix<double, 2 * dim, 1>::Zero(2 * info.A.rows());
        Eigen::Index counter = T - begin;
        integrate_times(stepper, closed_loop{ info }, zeros, T, end, 1e-5,
            [&](const state_t& state, double t) {
                ret.col(counter++) = state;
            });
    }
    Eigen::Matrix<double, 2 * dim, 2 * dim> boundary;
    boundary.topRows<dim>(info.A.rows()) = cache.LTI_STM[0].template topRows<dim>(info.A.rows());
    boundary.bottomRows<dim>(info.A.rows()) = cache.LTI_STM[T - 1 - begin].template bottomRows<dim>(info.A.rows());
    Eigen::Matrix<double, 2 * dim, 1> shift;
    shift.topRows<dim>(info.A.rows()) = x - ret.topLeftCorner(info.A.rows(), 1);
    shift.bottomRows<dim>(info.A.rows()).setZero();
    Eigen::Matrix<double, 2 * dim, 1> corr = boundary.colPivHouseholderQr().solve(shift);
    for (std::size_t t = 0; t != end - begin; ++t)
        ret.col(t) += cache.LTI_STM[t] * corr;
    return ret;
}
std::pair<Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>, double> solver::ode_init(const Eigen::Matrix<double, dim, 1>& x, std::size_t len) const {
    fixed_point_check check{ __func__ };
    auto [mu_short, err] = fixed_point<Eigen::Matrix<double, mu_dim, Eigen::Dynamic>>(
        [&](const mu_matrix& mu) -> mu_matrix {
            Eigen::Matrix<double, mu_dim, Eigen::Dynamic> rhs =
                mult * qint_LTI(x, mu, info.t.cbegin() + len, info.t.cbegin());
            return info.eta_inv(rhs);
        }, mu_matrix::Zero(mu_dim, len), 220, 1e-8, 1e-8 * len, check);
    return std::make_pair(qint_LTI(x, mu_short, info.t.cbegin() + len, info.t.cbegin(), info.t.cend()), err);
}
Eigen::Matrix<double, 2 * dim, 1>
solver::ode_rhs(const Eigen::Matrix<double, dim, 1>& x, Eigen::Matrix<double, 2 * dim, 1> state, double T) const {
    Eigen::Matrix<double, 2 * dim, 1 + dim> init;
    const Eigen::Index _dim = info.A.rows();
    init.col(0) = state;
    init.rightCols<dim>(_dim) = Eigen::Matrix<double, 2 * dim, 2 * dim>::Identity(_dim, _dim).rightCols<dim>(_dim);

    {
    typedef Eigen::Matrix<double, 2 * dim, 1> state_t;
    using namespace boost::numeric::odeint;
    auto stepper = make_controlled<runge_kutta_dopri5<state_t>>(1e-10, 1e-10);
    integrate_adaptive(stepper, closed_loop{ info }, state, info.t.front(), T, 1e-5);
    }
    // state is now at final time
    Eigen::Matrix<double, 2 * dim, 1> shift;
    closed_loop{ info }(state, shift, 0);
    {
    using namespace boost::numeric::odeint;
    auto stepper = make_controlled<runge_kutta_dopri5<Eigen::Matrix<double, 2 * dim, 1 + dim>>>(1e-10, 1e-10);
    integrate_adaptive(stepper, closed_loop_with_der{ info }, init, info.t.front(), T, 1e-5);
    }
    // Eigen::Matrix<double, 2 * dim, 2 * dim> boundary;
    // boundary.topRows<dim>(info.A.rows()).setIdentity(); // der(0).topRows<dim>();
    // boundary.bottomRows<dim>(info.A.rows()) = init.bottomRightCorner<dim, 2 * dim>(_dim, 2 * _dim);
    Eigen::Matrix<double, 2 * dim, 1> ret(_dim);
    ret.topRows<dim>(_dim).setZero();
    ret.bottomRows<dim>(_dim) = init.bottomRightCorner<dim, dim>(_dim, _dim).colPivHouseholderQr().solve(-shift.bottomRows<dim>(info.A.rows()));
    return ret;
}
Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>
solver::ode_next_euler(const Eigen::Matrix<double, dim, 1>& x, const Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>& z_Pzq, Eigen::Index T) const {
    Eigen::Matrix<double, 2 * dim, 1 + 2 * dim> init;
    const Eigen::Index _dim = info.A.rows();
    init.col(0) = z_Pzq.col(0);
    init.rightCols<2 * dim>(2 * _dim).setIdentity();
    Eigen::Matrix<double, 2 * dim, 1> shift;
    closed_loop{ info }(z_Pzq.col(T), shift, 0);
    using namespace boost::numeric::odeint;
    auto stepper = make_dense_output<runge_kutta_dopri5<Eigen::Matrix<double, 2 * dim, 1 + 2 * dim>>>(1e-10, 1e-10);
    Eigen::Index counter = 0;
    Eigen::Matrix<double, 2 * dim, Eigen::Dynamic> z_lambda(2 * dim, T + 2);
    std::vector<Eigen::Matrix<double, 2 * dim, 2 * dim>> der(T + 2);
    integrate_times(stepper, closed_loop_with_der{ info }, init, info.t.cbegin(), info.t.cbegin() + T + 2, 1e-5,
        [&](const auto& state, double) {
            z_lambda.col(counter) = state.col(0);
            der[counter++] = state.template rightCols<2 * dim>(2 * _dim);
        });
    Eigen::Matrix<double, 2 * dim, 2 * dim> boundary;
    boundary.topRows<dim>(info.A.rows()) = der[0].topRows<dim>(info.A.rows());
    boundary.bottomRows<dim>(info.A.rows()) = der[T].bottomRows<dim>(info.A.rows());
    shift.topRows<dim>(info.A.rows()).setZero();
    shift.bottomRows<dim>(info.A.rows()) = -shift.bottomRows<dim>(info.A.rows()).eval();
    Eigen::Matrix<double, 2 * dim, 1> corr = (info.t[T + 1] - info.t[T]) * boundary.colPivHouseholderQr().solve(shift);
    for (Eigen::Index t = 0; t <= T + 1; ++t)
        z_lambda.col(t) += der[t] * corr;
    return z_lambda;
}
void write_mu(const PDE_info& info, mu_matrix& mu, const Eigen::Matrix<double, dim + dim, Eigen::Dynamic>& pq, const RS_interp& RS) {
    const Eigen::Index N = pq.cols() - 1;
    for (Eigen::Index j = 0; j != N; ++j) {
        auto p0 = pq.col(j).topRows<dim>();
        auto q0 = pq.col(j).bottomRows<dim>();
        const auto& [R0, S0] = RS[2 * j];
        const auto& [R1, S1] = RS[2 * j + 1];
        Eigen::Matrix<double, mu_dim, 1> y;
        y.topRows<M_dim>() = info.M_0 * (p0 + R0 * q0);
        y.bottomRows<L_dim>() = info.L_0.transpose() * (S0 * q0);
        mu.col(2 * j) = info.eta_inv(y);
        Eigen::Matrix<double, dim, 1> p1 = 0.5 * (p0 + pq.col(j + 1).topRows<dim>());
        Eigen::Matrix<double, dim, 1> q1 = 0.5 * (q0 + pq.col(j + 1).bottomRows<dim>());
        y.topRows<M_dim>() = info.M_0 * (p1 + R1 * q1);
        y.bottomRows<L_dim>() = info.L_0.transpose() * (S1 * q1);
        mu.col(2 * j + 1) = info.eta_inv(y);
    }
    auto p0 = pq.col(N).topRows<dim>();
    auto q0 = pq.col(N).bottomRows<dim>();
    const auto& [R0, S0] = RS[2 * N];
    Eigen::Matrix<double, mu_dim, 1> y;
    y.topRows<M_dim>() = info.M_0 * (p0 + R0 * q0);
    y.bottomRows<L_dim>() = info.L_0.transpose() * (S0 * q0);
    mu.col(2 * N) = info.eta_inv(y);
}
std::tuple<double, Eigen::Matrix<double, dim, 1>, Eigen::Matrix<double, dim, 1>>
solver::solve(const Eigen::Matrix<double, dim, 1>& x, std::span<const double>, std::span<const std::size_t> N, std::span<const double> errs) const {
    std::vector<RS_interp> RS;
    RS.reserve(info.t.size() - 1);
    for (std::size_t i = 0; i != info.t.size() - 1; ++i)
        RS.emplace_back(info, info.t[i], info.t[i + 1], 2 * N[i]);
    using namespace boost::numeric::odeint;
    typedef Eigen::Matrix<double, 2 * dim, 1> state_t;
    auto stepper = runge_kutta4<state_t>();
    state_t tmp;
    tmp.topRows<dim>() = x;
    tmp.bottomRows<dim>().setZero();
    Eigen::Matrix<double, dim, Eigen::Dynamic> P(dim, info.t.size());
    {
    int i = 0;
    integrate_n_steps(stepper, closed_loop{info}, tmp, info.t.front(), info.t[1] - info.t[0], info.t.size() - 1, [&](const state_t& state, double) {
        P.col(i++) = state.topRows<dim>();
    });
    }
    Eigen::Matrix<double, dim, Eigen::Dynamic> Q(dim, info.t.size());
    Q.setZero(); // FIXME
    std::vector<mu_matrix> mu;
    for (int i = 0; i != info.t.size() - 1; ++i)
        mu.push_back(mu_matrix::Zero(mu_dim, 2 * N[i] + 1));
    bool down = false;
    for (double err : errs) {
        int k = info.t.size() - 1;
        for (std::size_t counter = 0; counter < 100 * info.t.size() * info.t.size(); ++counter) {
            if (down) {
                down = false;
                ++k;
            }
            else if (k > 0)
                --k;
            else
                break;
            RS_backward pq_ode{info, RS[k], mu[k]};
            Eigen::Matrix<double, dim + dim, Eigen::Dynamic> pq(dim + dim, N[k] + 1);
            for (uint_fast8_t i = 0; i != 255; ++i) {
                mu_matrix old_mu = mu[k];
                tmp.topRows<dim>() = P.col(k);
                tmp.bottomRows<dim>().setZero();
                {
                int j = 0;
                integrate_n_steps(stepper, pq_ode, tmp, info.t[k], (info.t[k + 1] - info.t[k]) / N[k], N[k], [&](const state_t& state, double) {
                    pq.col(j++) = state;
                });
                }
                auto p = pq.col(N[k]).topRows<dim>();
                if (k + 1 == info.t.size() - 1)
                    Q.col(k + 1).setZero(); // FIXME
                else if ((p + RS[k].R_end() * Q.col(k + 1) - P.col(k + 1)).norm() > err) {
                    P.col(k + 1) = p + RS[k].R_end() * Q.col(k + 1);
                    down = true;
                    break;
                }
                pq.bottomRows<dim>().colwise() += (Q.col(k + 1) - pq.col(N[k]).bottomRows<dim>()).eval();
                write_mu(info, mu[k], pq, RS[k]);
                Q.col(k) = RS[k].S_begin() * pq.col(0).bottomRows<dim>();
                double diff = (mu[k] - old_mu).norm();
                if (diff < err)
                    break;
            }
        }
    }
    Eigen::Matrix<double, 2 * dim + 1, 1> state;
    state.topRows<dim>() = x;
    state.middleRows<dim>(dim) = Q.col(0);
    state[2 * dim] = 0;
    {
    auto stepper = make_controlled<runge_kutta_dopri5<decltype(state)>>(1e-10, 1e-10);
    integrate_adaptive(stepper, closed_loop_with_r{ info },
        state, info.t.front(), info.t.back(), 1e-5);
    }
    return std::make_tuple(state[2 * dim], Q.col(0), state.middleRows<dim>(dim).eval());
}
Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>
solver::solve_times(const Eigen::Matrix<double, dim, 1>& x, std::size_t T) const {
    Eigen::Matrix<double, 2 * dim, Eigen::Dynamic> ret(2 * dim, info.t.size());
    auto [z_Pzq, err] = ode_init(x, T);
    ret.leftCols(T).colwise() = z_Pzq.col(0);
#if 0
    for (T; T != info.t.size(); ++T) {
        z_Pzq = ode_next_euler(x, z_Pzq, T - 1);
        ret.col(T) = z_Pzq.col(0);
    }
#else
    typedef Eigen::Matrix<double, 2 * dim, 1> state_t;
    using namespace boost::numeric::odeint;
    state_t state = z_Pzq.col(0);
    auto stepper = adams_bashforth_moulton<3, state_t>();
    Eigen::Index counter = T - 1;
    integrate_times(stepper, [&](const state_t& in, state_t& out, double t){ out = ode_rhs(x, in, t); }, state, info.t.cbegin() + T - 1, info.t.cend(), info.t.back() - info.t.front(), [&](const state_t& state, double t) {
        ret.col(counter++) = state;
    });
#endif
    return ret;
}

