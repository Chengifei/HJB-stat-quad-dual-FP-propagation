#include "solver.h"

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
std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::Matrix<double, dim, Eigen::Dynamic>>
solver::W(const Eigen::Matrix<double, dim, 1>& x, std::size_t T) const {
    const Eigen::Index _dim = info.A.rows();
    auto inits = solve_times(x, T);
    Eigen::VectorXd ret(info.t.size());
    Eigen::VectorXd err(info.t.size());
    typedef Eigen::Matrix<double, 2 * dim + 1, 1> state_t;
    using namespace boost::numeric::odeint;
    auto stepper = make_dense_output<runge_kutta_dopri5<state_t>>(1e-10, 1e-10);
    for (std::size_t i = 0; i != info.t.size(); ++i) {
        state_t state;
        state.template topRows<dim + dim>(_dim + _dim) = inits.col(i);
        state[_dim + _dim] = 0;
        err[i] = (state.topRows<dim>(_dim) - x).norm();
        integrate_adaptive(stepper, closed_loop_with_r{ info },
            state, info.t.front(), info.t[i], 1e-5);
        err[i] += (state.middleRows<dim>(_dim) - Eigen::Matrix<double, dim, 1>::Zero(_dim, 1)).norm();
        ret[i] = state[_dim + _dim];
    }
    return std::make_tuple(std::move(ret), std::move(err), inits.bottomRows<dim>(_dim));
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
std::pair<Eigen::Matrix<double, 2 * dim, Eigen::Dynamic>, Eigen::VectorXd>
solver::solve(const Eigen::Matrix<double, dim, 1>& x, std::size_t T) const {
    auto [z_Pzq, err] = ode_init(x, T);
    Eigen::VectorXd res(info.t.size());
    std::fill(res.begin(), res.begin() + T, err);
    for (T; T != info.t.size(); ++T) {
        z_Pzq = ode_next_euler(x, z_Pzq, T - 1);
        mu_matrix mu = info.eta_inv(mult * z_Pzq.leftCols(T + 1));
        res[T] = (mu - info.eta_inv(mult * qint_LTI(x, mu, info.t.cbegin() + T + 1, info.t.cbegin()))).norm();
    }
    return std::make_pair(z_Pzq.leftCols(T), std::move(res));
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

