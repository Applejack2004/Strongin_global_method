#include <cmath>
#include <vector>
#include <limits>
#include <fstream>
#include <chrono>
#include <iostream>
#include <random>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iomanip>
#include <string>
#include <filesystem>
#include <clocale>
#include <sstream>
#include <map>
#include <functional>
#include <memory>

// Подключение всех необходимых классов задач
#include "HillProblem.hpp"
#include "ShekelProblem.hpp"
#include "HillProblemFamily.hpp"
#include "ShekelProblemFamily.hpp"
#include "grishagin_function.hpp"
#include "GrishaginProblemFamily.hpp"
#include "OptSqConstrProblem.hpp"
#include "HillConstrainedProblem.hpp"
#include "ShekelConstrainedProblem.hpp"
#include "GKLSConstrainedProblem.hpp"
#include "GKLSProblem.hpp"
#include <cstdlib>

// Подключение библиотеки DIRECT
#include "direct.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

struct LocalizationSettings {
    bool enabled = false;
    int delay = 0;
    double alpha = 10.0;
    int global_ratio = 1;
    int local_ratio = 1;
};

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class IConstrainedProblem;

IConstrainedOptProblem* g_current_constrained_problem = nullptr;
IOptProblem* g_current_unconstrained_problem = nullptr;
IConstrainedProblem* g_current_generic_constrained_problem = nullptr;
double g_penalty_coefficient = 1000.0;

std::ofstream g_direct_points_file;
std::vector<double> g_direct_lb;
std::vector<double> g_direct_ub;
bool g_save_direct_points = false;

int g_direct_force_stop = 0;
double g_direct_eps = 0.0;
std::vector<double> g_direct_optimum;
std::vector<double> g_direct_stopped_x;
double g_direct_stopped_f = 0.0;

// ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ ДЛЯ РЕЖИМА СТРОГОГО БЮДЖЕТА (ДЕБАГ)
std::ofstream g_debug_convergence_file;
double g_current_best_val = std::numeric_limits<double>::infinity();
int g_eval_counter = 0;
bool g_strict_budget_mode = false;

void setup_direct_stopping(const std::vector<double>& lb, const std::vector<double>& ub, double eps, const std::vector<double>& opt) {
    g_direct_lb = lb;
    g_direct_ub = ub;
    g_direct_eps = eps;
    g_direct_optimum = opt;
    g_direct_force_stop = 0;
    g_direct_stopped_x.clear();
    g_direct_stopped_f = 0.0;
}

void open_direct_log(const std::string& filename) {
    g_direct_points_file.open(filename);
    g_save_direct_points = true;
}

void close_direct_log() {
    if (g_direct_points_file.is_open()) {
        g_direct_points_file.close();
    }
    g_save_direct_points = false;
}

bool check_direct_stop(const std::vector<double>& point) {
    if (g_direct_optimum.empty() || g_direct_lb.empty() || g_direct_ub.empty()) {
        return false;
    }
    bool all_close = true;
    double effective_eps = g_direct_eps / 2.0;

    for (size_t i = 0; i < point.size(); ++i) {
        double range = g_direct_ub[i] - g_direct_lb[i];
        if (range > 1e-12) {
            double diff_normalized = std::abs(point[i] - g_direct_optimum[i]) / range;
            if (diff_normalized > effective_eps) {
                all_close = false;
                break;
            }
        }
    }
    return all_close;
}

double penalty_function_for_direct(int n, const double* x, int* undefined_flag, void* data);
double generic_penalty_function_for_direct(int n, const double* x, int* undefined_flag, void* data);

double objective_function_for_direct(int n, const double* x, int* undefined_flag, void* data)
{
    (void)n;
    (void)undefined_flag;

    if (!g_current_unconstrained_problem) return 0.0;

    vector<double> point(x, x + n);
    double final_z = g_current_unconstrained_problem->ComputeFunction(point);

    if (g_strict_budget_mode) {
        g_eval_counter++;
        if (final_z < g_current_best_val) {
            g_current_best_val = final_z;
        }
        if (g_debug_convergence_file.is_open()) {
            g_debug_convergence_file << "итерация " << g_eval_counter << ": " << std::fixed << std::setprecision(8) << g_current_best_val << "\n";
        }
    } else {
        if (check_direct_stop(point)) {
            g_direct_force_stop = 1;
            g_direct_stopped_x = point;
            g_direct_stopped_f = final_z;
        }
    }

    if (g_save_direct_points && g_direct_points_file.is_open()) {
        for (int d = 0; d < n; ++d) {
            double mapped = point[d];
            if (g_direct_ub.size() == (size_t)n && (g_direct_ub[d] - g_direct_lb[d]) > 1e-12) {
                mapped = (point[d] - g_direct_lb[d]) / (g_direct_ub[d] - g_direct_lb[d]);
            }
            g_direct_points_file << mapped << " ";
        }
        g_direct_points_file << final_z << " 1\n";
    }
    return final_z;
}

int n1, nexp, l, iq, iu[10], iv[10];

void node(int is_param) {
    int current_n_dim, i, j_idx, k1, k2, iff_val;
    current_n_dim = n1 + 1;
    if (is_param == 0) {
        l = n1;
        for (i = 0; i < current_n_dim; i++) {
            iu[i] = -1;
            iv[i] = -1;
        }
    }
    else if (is_param == (nexp - 1)) {
        l = n1;
        iu[0] = 1;
        iv[0] = 1;
        for (i = 1; i < current_n_dim; i++) {
            iu[i] = -1;
            iv[i] = -1;
        }
        iv[n1] = 1;
    }
    else {
        iff_val = nexp;
        k1 = -1;
        for (i = 0; i < current_n_dim; i++) {
            iff_val = iff_val / 2;
            if (is_param >= iff_val) {
                if ((is_param == iff_val) && (is_param != 1)) {
                    l = i;
                    iq = -1;
                }
                is_param = is_param - iff_val;
                k2 = 1;
            }
            else {
                k2 = -1;
                if ((is_param == (iff_val - 1)) && (is_param != 0)) {
                    l = i;
                    iq = 1;
                }
            }
            j_idx = -k1 * k2;
            iv[i] = j_idx;
            iu[i] = j_idx;
            k1 = k2;
        }
        iv[l] = iv[l] * iq;
        iv[n1] = -iv[n1];
    }
}

void mapd(double x_input, int m_order, double* y_output, int n_dim, int key_type) {
    double mne_val, dd_val, dr_val;
    double p_calc, r_calc;
    int iw_arr[11];
    int it_val, is_val, i_loop, j_loop, k_idx = 0;

    p_calc = 0.0;
    n1 = n_dim - 1;
    if (n_dim <= 0 || n_dim >= 10) {
        for (i_loop = 0; i_loop < n_dim; i_loop++) {
            y_output[i_loop] = 0.5;
        }
        return;
    }

    for (nexp = 1, i_loop = 0; i_loop < n_dim; nexp *= 2, i_loop++) {}

    if (nexp <= 0 && n_dim > 0) nexp = 1;

    double d_for_key_logic = x_input;
    r_calc = 0.5;
    it_val = 0;
    dr_val = static_cast<double>(nexp);
    for (mne_val = 1.0, i_loop = 0; i_loop < m_order; mne_val *= dr_val, i_loop++) {}

    for (i_loop = 0; i_loop < n_dim; i_loop++) {
        iw_arr[i_loop] = 1;
        y_output[i_loop] = 0.0;
    }

    if (key_type == 2) {
        if (std::abs(mne_val) > 1e-12) {
            d_for_key_logic = d_for_key_logic * (1.0 - 1.0 / mne_val);
        } else {
            d_for_key_logic = (x_input == 1.0) ? 1.0 : 0.0;
        }
        k_idx = 0;
    }
    else if (key_type > 2) {
        double temp_d_for_key_gt2 = x_input;
        dr_val = mne_val / nexp;
        dr_val = dr_val - fmod(dr_val, 1.0);
        dd_val = mne_val - dr_val;
        dr_val = temp_d_for_key_gt2 * dd_val;
        dd_val = dr_val - fmod(dr_val, 1.0);

        if (nexp - 1 != 0) {
            dr_val = dd_val + (dd_val - 1.0) / (nexp - 1.0);
        } else {
            dr_val = dd_val;
        }

        dd_val = dr_val - fmod(dr_val, 1.0);

        if (std::abs(mne_val) > 1e-12) {
            d_for_key_logic = dd_val * (1.0 / mne_val);
        } else {
            d_for_key_logic = 0.0;
        }
    }

    double d_iter_for_digits = (key_type == 1) ? x_input : d_for_key_logic;
    double final_d_remainder = 0.0;

    for (j_loop = 0; j_loop < m_order; j_loop++) {
        iq = 0;
        if (std::abs(x_input - 1.0) < 1e-9 && key_type != 2) {
            is_val = nexp - 1;
            if (j_loop == m_order - 1) {
                final_d_remainder = 0.0;
            }
        }
        else {
            d_iter_for_digits = d_iter_for_digits * static_cast<double>(nexp);
            is_val = static_cast<int>(floor(d_iter_for_digits));
            if (is_val >= nexp) is_val = nexp - 1;
            if (is_val < 0) is_val = 0;

            d_iter_for_digits = d_iter_for_digits - static_cast<double>(is_val);
            if (j_loop == m_order - 1) final_d_remainder = d_iter_for_digits;
        }

        node(is_val);

        int temp_i;
        temp_i = iu[0]; iu[0] = iu[it_val]; iu[it_val] = temp_i;
        temp_i = iv[0]; iv[0] = iv[it_val]; iv[it_val] = temp_i;

        if (l == 0 && it_val != 0) l = it_val;
        else if (l == it_val) l = 0;

        if ((iq > 0) || ((iq == 0) && (is_val == 0))) k_idx = l;
        else if (iq < 0) k_idx = (it_val == n1) ? 0 : n1;

        r_calc = r_calc * 0.5;
        it_val = l;
        for (i_loop = 0; i_loop < n_dim; i_loop++) {
            iu[i_loop] = iu[i_loop] * iw_arr[i_loop];
            iw_arr[i_loop] = -iv[i_loop] * iw_arr[i_loop];
            p_calc = r_calc * static_cast<double>(iu[i_loop]);
            p_calc = p_calc + y_output[i_loop];
            y_output[i_loop] = p_calc;
        }
    }

    if (key_type == 2) {
        int sign_val = (is_val == (nexp - 1)) ? -1 : 1;
        if (k_idx < 0 || k_idx >= 10) k_idx = 0;
        p_calc = 2.0 * static_cast<double>(sign_val) * static_cast<double>(iu[k_idx]) * r_calc * final_d_remainder;
        y_output[k_idx] = y_output[k_idx] - p_calc;
    }
    else if (key_type == 3) {
        for (i_loop = 0; i_loop < n_dim; i_loop++) {
            if (i_loop >= 10) continue;
            p_calc = r_calc * static_cast<double>(iu[i_loop]);
            p_calc = p_calc + y_output[i_loop];
            y_output[i_loop] = p_calc;
        }
    }
}

void resetMappingGlobals(int) {
    n1 = 0; nexp = 1; l = 0; iq = 0;
    for (int i = 0; i < 10; ++i) { iu[i] = 0; iv[i] = 0; }
}

std::vector<double> peanoMapping(double x_param, int m_order, int n_dim, int key_type) {
    if (n_dim <= 0 || n_dim >= ((sizeof(iu) / sizeof(iu[0])) - 1)) {
        return std::vector<double>(n_dim > 0 ? n_dim : 1, 0.5);
    }
    resetMappingGlobals(n_dim);
    std::vector<double> y_result(n_dim);
    mapd(x_param, m_order, y_result.data(), n_dim, key_type);
    for (int i = 0; i < n_dim; ++i) y_result[i] += 0.5;
    return y_result;
}

struct Point {
    double x_param;
    double z_value;
    std::vector<double> y_coords;
};

class FunctionInterface {
public:
    virtual double ComputeFunction(const std::vector<double>& x) const = 0;
    virtual std::vector<double> GetOptimumPoint() const = 0;
    virtual double GetOptimumValue() const = 0;
    virtual ~FunctionInterface() {}
};

class IConstrainedProblem {
public:
    virtual double ComputeObjective(const std::vector<double>& x) const = 0;
    virtual double ComputeConstraint(int constraint_idx, const std::vector<double>& x) const = 0;
    virtual int GetConstraintsCount() const = 0;
    virtual std::vector<double> GetOptimumPoint() const { return {}; }
    virtual ~IConstrainedProblem() {}
};

class Generic1DProblem : public IConstrainedProblem {
private:
    std::function<double(double)> objective_func;
    std::vector<std::function<double(double)>> constraint_funcs;

public:
    Generic1DProblem(std::function<double(double)> obj_func) : objective_func(obj_func) {}
    void AddConstraint(std::function<double(double)> constr_func) { constraint_funcs.push_back(constr_func); }
    double ComputeObjective(const std::vector<double>& x) const override { return objective_func(x[0]); }
    double ComputeConstraint(int constraint_idx, const std::vector<double>& x) const override {
        return constraint_funcs[constraint_idx](x[0]);
    }
    int GetConstraintsCount() const override { return static_cast<int>(constraint_funcs.size()); }
};

class Example51Problem : public IConstrainedProblem {
public:
    double ComputeObjective(const std::vector<double>& u) const override {
        double y1 = u[0] * 4.0;
        double y2 = -1.0 + u[1] * 4.0;
        double term1 = -1.5 * y1 * y1 * exp(1.0 - y1 * y1 - 20.25 * pow(y1 - y2, 2));
        double term2 = -pow(0.5 * (y1 - 1.0) * (y2 - 1.0), 4.0) * exp(2.0 - pow(0.5 * (y1 - 1.0), 4.0) - pow(y2 - 1.0, 4.0));
        return term1 + term2;
    }
    double ComputeConstraint(int idx, const std::vector<double>& u) const override {
        double y1 = u[0] * 4.0;
        double y2 = -1.0 + u[1] * 4.0;
        if (idx == 0) return 0.01 * (pow(y1 - 2.2, 2) + pow(y2 - 1.2, 2) - 2.25);
        else if (idx == 1) return 100.0 * (1.0 - pow(y1 - 2.0, 2) / 1.44 - pow(0.5 * y2, 2));
        else if (idx == 2) return 10.0 * (y2 - 1.5 - 1.5 * sin(6.283 * (y1 - 1.75)));
        return 0.0;
    }
    int GetConstraintsCount() const override { return 3; }
};

double penalty_function_for_direct(int n, const double* x, int* undefined_flag, void* data)
{
    (void)undefined_flag; (void)data;
    if (!g_current_constrained_problem) return 0.0;
    vector<double> point(x, x + n);
    double objective_value = g_current_constrained_problem->ComputeFunction(point);
    double penalty = 0.0;
    int num_constraints = g_current_constrained_problem->GetConstraintsNumber();
    int v = num_constraints + 1;

    for (int i = 0; i < num_constraints; ++i) {
        double constraint_value = g_current_constrained_problem->ComputeConstraint(i, point);
        if (constraint_value > 0) {
            if (v == num_constraints + 1) v = i + 1;
            penalty += constraint_value * constraint_value;
        }
    }
    double final_z = objective_value + g_penalty_coefficient * penalty;

    if (g_strict_budget_mode) {
        g_eval_counter++;
        if (penalty == 0.0 && objective_value < g_current_best_val) {
            g_current_best_val = objective_value;
        }
        if (g_debug_convergence_file.is_open()) {
            if (g_current_best_val == std::numeric_limits<double>::infinity()) {
                g_debug_convergence_file << "итерация " << g_eval_counter << ": N/A\n";
            } else {
                g_debug_convergence_file << "итерация " << g_eval_counter << ": " << std::fixed << std::setprecision(8) << g_current_best_val << "\n";
            }
        }
    } else {
        if (penalty == 0.0 && check_direct_stop(point)) {
            g_direct_force_stop = 1; g_direct_stopped_x = point; g_direct_stopped_f = final_z;
        }
    }

    if (g_save_direct_points && g_direct_points_file.is_open()) {
        for (int d = 0; d < n; ++d) {
            double mapped = point[d];
            // ИСПРАВЛЕНИЕ: Масштабируем до [0,1] ТОЛЬКО если размерность больше 1
            if (n > 1 && g_direct_ub.size() == (size_t)n && (g_direct_ub[d] - g_direct_lb[d]) > 1e-12) {
                mapped = (point[d] - g_direct_lb[d]) / (g_direct_ub[d] - g_direct_lb[d]);
            }
            g_direct_points_file << mapped << " ";
        }
        g_direct_points_file << final_z << " " << v << "\n";
    }
    return final_z;
}

double generic_penalty_function_for_direct(int n, const double* x, int* undefined_flag, void* data)
{
    (void)undefined_flag; (void)data;
    if (!g_current_generic_constrained_problem) return 0.0;
    vector<double> point(x, x + n);
    double objective_value = g_current_generic_constrained_problem->ComputeObjective(point);
    double penalty = 0.0;
    int num_constraints = g_current_generic_constrained_problem->GetConstraintsCount();
    int v = num_constraints + 1;

    for (int i = 0; i < num_constraints; ++i) {
        double constraint_value = g_current_generic_constrained_problem->ComputeConstraint(i, point);
        if (constraint_value > 0) {
            if (v == num_constraints + 1) v = i + 1;
            penalty += constraint_value * constraint_value;
        }
    }
    double final_z = objective_value + g_penalty_coefficient * penalty;

    if (g_strict_budget_mode) {
        g_eval_counter++;
        if (penalty == 0.0 && objective_value < g_current_best_val) {
            g_current_best_val = objective_value;
        }
        if (g_debug_convergence_file.is_open()) {
            if (g_current_best_val == std::numeric_limits<double>::infinity()) {
                g_debug_convergence_file << "итерация " << g_eval_counter << ": N/A\n";
            } else {
                g_debug_convergence_file << "итерация " << g_eval_counter << ": " << std::fixed << std::setprecision(8) << g_current_best_val << "\n";
            }
        }
    } else {
        if (penalty == 0.0 && check_direct_stop(point)) {
            g_direct_force_stop = 1; g_direct_stopped_x = point; g_direct_stopped_f = final_z;
        }
    }

    if (g_save_direct_points && g_direct_points_file.is_open()) {
        for (int d = 0; d < n; ++d) {
            double mapped = point[d];
            // ИСПРАВЛЕНИЕ: Масштабируем до [0,1] ТОЛЬКО если размерность больше 1
            if (n > 1 && g_direct_ub.size() == (size_t)n && (g_direct_ub[d] - g_direct_lb[d]) > 1e-12) {
                mapped = (point[d] - g_direct_lb[d]) / (g_direct_ub[d] - g_direct_lb[d]);
            }
            g_direct_points_file << mapped << " ";
        }
        g_direct_points_file << final_z << " " << v << "\n";
    }
    return final_z;
}

// === ФАБРИКА ЗАДАЧ ===
std::unique_ptr<IConstrainedProblem> ProblemFactory(const std::string& problem_name) {
    if (problem_name == "Example3.2") {
        auto objective =[](double x) { return cos(18.0 * x - 3.0) * sin(10.0 * x - 7.0) + 1.0; };
        auto problem = std::make_unique<Generic1DProblem>(objective);
        problem->AddConstraint([](double x) { return exp(-x / 2.0) * sin(6.0 * x - 1.5); });
        problem->AddConstraint([](double x) { return sin(4.0 * x - 2.2) + cos(6.0 * x - 2.9); });
        problem->AddConstraint([](double x) { return std::abs(x) * sin(2.0 * M_PI * x - 0.5); });
        return problem;
    }
    else if (problem_name == "OptSqConstr1D") {
        auto problem_impl = std::make_shared<TOptSqConstrProblem>(
            1, vector<double>{-1.0}, vector<double>{1.0}, vector<double>{0.0}, 0.0, cptNormal, 1.0, 1
        );
        auto wrapper = std::make_unique<Generic1DProblem>([problem_impl](double x) {
            return problem_impl->ComputeFunction({ x });
        });
        for (int i = 0; i < problem_impl->GetConstraintsNumber(); ++i) {
            wrapper->AddConstraint([problem_impl, i](double x) { return problem_impl->ComputeConstraint(i, { x }); });
        }
        return wrapper;
    }
    return nullptr;
}
// ===================================

template <typename T_ProblemType> class Minimizer {
private:
    std::vector<double> searchLeftBound_param;
    std::vector<double> searchRightBound_param;
    double epsilon_val;
    double r_strongin_parameter;
    const T_ProblemType& function_instance;
    std::ostream& logger;

    int max_iterations_limit_member;
    int iteration_counter;
    int exit_main_criteria_count;
    int exit_test_criteria_count;

    std::vector<Point> trial_points_list;

    bool use_peano_mapping_flag;
    int peano_mapping_m_order;
    int problem_actual_dimension;
    int peano_mapping_key_type;

    double current_r_mu_estimate;
    double current_max_mu_estimate;
    LocalizationSettings loc_settings;
    double inv_dim;

    double getHolderDelta(double x_left, double x_right) const {
        double dx_param_abs = std::abs(x_right - x_left);
        if (!use_peano_mapping_flag || problem_actual_dimension <= 0) return dx_param_abs;
        if (problem_actual_dimension == 1 && use_peano_mapping_flag) return dx_param_abs;
        if (dx_param_abs < 1e-15) return 1e-15;
        return std::pow(dx_param_abs, inv_dim);
    }

    void performFirstIteration() {
        Point p_left;
        p_left.x_param = searchLeftBound_param[0];
        if (!use_peano_mapping_flag) p_left.y_coords = { p_left.x_param };
        else p_left.y_coords = peanoMapping(p_left.x_param, peano_mapping_m_order, problem_actual_dimension, peano_mapping_key_type);
        p_left.z_value = function_instance.ComputeFunction(p_left.y_coords);
        trial_points_list.push_back(p_left);

        Point p_right;
        p_right.x_param = searchRightBound_param[0];
        if (!use_peano_mapping_flag) p_right.y_coords = { p_right.x_param };
        else p_right.y_coords = peanoMapping(p_right.x_param, peano_mapping_m_order, problem_actual_dimension, peano_mapping_key_type);
        p_right.z_value = function_instance.ComputeFunction(p_right.y_coords);
        trial_points_list.push_back(p_right);

        if (use_peano_mapping_flag && (searchLeftBound_param[0] <= 0.5 && searchRightBound_param[0] >= 0.5) && problem_actual_dimension > 0) {
            Point center_point;
            center_point.x_param = 0.5;
            center_point.y_coords = peanoMapping(center_point.x_param, peano_mapping_m_order, problem_actual_dimension, peano_mapping_key_type);
            center_point.z_value = function_instance.ComputeFunction(center_point.y_coords);
            trial_points_list.push_back(center_point);
        }

        std::sort(trial_points_list.begin(), trial_points_list.end(),[](const Point& a, const Point& b) {
            return a.x_param < b.x_param;
        });
    }

    double computeIntervalCharacteristicR(const Point& p1, const Point& p2) const {
        double delta_i = getHolderDelta(p1.x_param, p2.x_param);
        if (delta_i < 1e-12) return -std::numeric_limits<double>::infinity();
        double dz = p2.z_value - p1.z_value;
        double reliable_r_mu = std::max(current_r_mu_estimate, r_strongin_parameter);
        return reliable_r_mu * delta_i + (dz * dz) / (reliable_r_mu * delta_i) - 2.0 * (p1.z_value + p2.z_value);
    }

    Point generateNewTrialPoint(size_t interval_index_left) {
        const Point& p_left = trial_points_list[interval_index_left];
        const Point& p_right = trial_points_list[interval_index_left + 1];
        Point new_p;

        double xL = p_left.x_param;
        double xR = p_right.x_param;
        double zL = p_left.z_value;
        double zR = p_right.z_value;

        if (use_peano_mapping_flag && problem_actual_dimension > 0) {
            double mu_v_calc = std::max(current_max_mu_estimate, 1e-9);
            if (current_max_mu_estimate < 1e-9) mu_v_calc = 1.0;
            double term_dz_mu = std::abs(zR - zL) / mu_v_calc;
            double power_term = std::pow(term_dz_mu, static_cast<double>(problem_actual_dimension));

            int sign_dz = 0;
            if (zR - zL > 1e-9) sign_dz = 1;
            else if (zR - zL < -1e-9) sign_dz = -1;

            new_p.x_param = 0.5 * (xL + xR) - static_cast<double>(sign_dz) * (1.0 / (2.0 * r_strongin_parameter)) * power_term;
        }
        else {
            double reliable_r_mu_for_new_point = std::max(current_r_mu_estimate, r_strongin_parameter);
            new_p.x_param = 0.5 * (xL + xR) - (zR - zL) / (2.0 * reliable_r_mu_for_new_point);
        }

        double interval_width = xR - xL;
        double min_rel_step = 0.0001;
        double abs_min_step = 1e-10;
        double step_from_boundary = std::max(abs_min_step, interval_width * min_rel_step);

        if (interval_width <= 1e-12) new_p.x_param = xL + interval_width / 2.0;
        else new_p.x_param = std::max(xL + step_from_boundary, std::min(xR - step_from_boundary, new_p.x_param));

        if (!use_peano_mapping_flag) new_p.y_coords = { new_p.x_param };
        else new_p.y_coords = peanoMapping(new_p.x_param, peano_mapping_m_order, problem_actual_dimension, peano_mapping_key_type);

        new_p.z_value = function_instance.ComputeFunction(new_p.y_coords);
        return new_p;
    }

    bool checkStoppingConditions(const Point& p_left_of_interval, const Point& p_right_of_interval) {
        if (use_peano_mapping_flag && !trial_points_list.empty() && !g_strict_budget_mode) {
            auto min_trial_it = std::min_element(trial_points_list.begin(), trial_points_list.end(),[](const Point& t1, const Point& t2) {
                return t1.z_value < t2.z_value;
            });
            const Point& y_current_best_point = *min_trial_it;

            std::vector<double> known_optimum_y_coords;
            try { known_optimum_y_coords = function_instance.GetOptimumPoint(); } catch (...) {}

            if (!known_optimum_y_coords.empty() && known_optimum_y_coords.size() == static_cast<size_t>(problem_actual_dimension)) {
                bool all_coords_close_enough = true;
                for (int j = 0; j < problem_actual_dimension; ++j) {
                    if (std::abs(y_current_best_point.y_coords[j] - known_optimum_y_coords[j]) > epsilon_val) {
                        all_coords_close_enough = false;
                        break;
                    }
                }
                if (all_coords_close_enough) {
                    exit_test_criteria_count++;
                    return true;
                }
            }
        }

        if (!g_strict_budget_mode) {
            double holder_interval_length = getHolderDelta(p_left_of_interval.x_param, p_right_of_interval.x_param);
            if (holder_interval_length <= epsilon_val) {
                exit_main_criteria_count++;
                return true;
            }
        }
        return false;
    }

public:
    Minimizer(const std::vector<double>& problem_domain_a, const std::vector<double>& problem_domain_b, double eps,
        double r_val, const T_ProblemType& func, std::ostream& log_stream, int max_iter, LocalizationSettings loc = {})
        : epsilon_val(eps), r_strongin_parameter(r_val), function_instance(func), logger(log_stream),
        max_iterations_limit_member(max_iter), iteration_counter(0), exit_main_criteria_count(0),
        exit_test_criteria_count(0), use_peano_mapping_flag(false), problem_actual_dimension(1),
        current_r_mu_estimate(r_val), current_max_mu_estimate(1.0), loc_settings(loc) {

        searchLeftBound_param = problem_domain_a;
        searchRightBound_param = problem_domain_b;
        inv_dim = problem_actual_dimension > 0 ? 1.0 / static_cast<double>(problem_actual_dimension) : 1.0;
    }

    Minimizer(double peano_search_param_a, double peano_search_param_b, double eps, double r_val,
        const T_ProblemType& func, int mapping_m_order, int original_problem_dimension,
        int mapping_key, std::ostream& log_stream, int max_iter, LocalizationSettings loc = {})
        : epsilon_val(eps), r_strongin_parameter(r_val), function_instance(func), logger(log_stream),
        max_iterations_limit_member(max_iter), iteration_counter(0), exit_main_criteria_count(0),
        exit_test_criteria_count(0), use_peano_mapping_flag(true), peano_mapping_m_order(mapping_m_order),
        problem_actual_dimension(original_problem_dimension), peano_mapping_key_type(mapping_key),
        current_r_mu_estimate(r_val), current_max_mu_estimate(1.0), loc_settings(loc) {

        searchLeftBound_param = { peano_search_param_a };
        searchRightBound_param = { peano_search_param_b };
        inv_dim = problem_actual_dimension > 0 ? 1.0 / static_cast<double>(problem_actual_dimension) : 1.0;
    }

    ~Minimizer() {}

    std::vector<double> findMinimum() {
        trial_points_list.clear();
        trial_points_list.reserve(max_iterations_limit_member + 3);
        iteration_counter = 0;
        exit_main_criteria_count = 0;
        exit_test_criteria_count = 0;

        performFirstIteration();

        if (trial_points_list.size() < 2) {
            if (!trial_points_list.empty()) return trial_points_list[0].y_coords;
            return {};
        }

        for (int current_iter_loop = 0; current_iter_loop < this->max_iterations_limit_member; ++current_iter_loop) {
            iteration_counter = current_iter_loop + 1;
            current_max_mu_estimate = 0.0;
            for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
                double holder_delta_for_slope = getHolderDelta(trial_points_list[i].x_param, trial_points_list[i + 1].x_param);
                if (holder_delta_for_slope > 1e-12) {
                    double mu_val = std::abs(trial_points_list[i + 1].z_value - trial_points_list[i].z_value) / holder_delta_for_slope;
                    if (mu_val > current_max_mu_estimate) current_max_mu_estimate = mu_val;
                }
            }
            if (current_max_mu_estimate < 1e-9) current_max_mu_estimate = 1.0;
            current_r_mu_estimate = r_strongin_parameter * current_max_mu_estimate;

            bool is_local_step = false;
            if (loc_settings.enabled && iteration_counter > loc_settings.delay) {
                int cycle = loc_settings.global_ratio + loc_settings.local_ratio;
                int step_in_cycle = (iteration_counter - loc_settings.delay - 1) % cycle;
                if (step_in_cycle >= loc_settings.global_ratio) is_local_step = true;
            }

            double z_star = std::min_element(trial_points_list.begin(), trial_points_list.end(),[](const Point& a, const Point& b) {
                return a.z_value < b.z_value;
            })->z_value;

            double max_R_val = -std::numeric_limits<double>::infinity();
            size_t interval_to_split_idx = 0;

            for (size_t i = 0; i < trial_points_list.size() - 1; ++i) {
                double R = computeIntervalCharacteristicR(trial_points_list[i], trial_points_list[i + 1]);

                if (is_local_step && R != -std::numeric_limits<double>::infinity()) {
                    double mu = std::max(current_max_mu_estimate, 1.0);
                    double zL = trial_points_list[i].z_value;
                    double zR = trial_points_list[i+1].z_value;
                    double delta_loc = std::sqrt(std::max(zL - z_star, 0.0) * std::max(zR - z_star, 0.0)) / mu;
                    R *= 1.0 / (delta_loc + std::pow(1.5, -loc_settings.alpha));
                }

                if (R > max_R_val) {
                    max_R_val = R;
                    interval_to_split_idx = i;
                }
            }

            if (checkStoppingConditions(trial_points_list[interval_to_split_idx], trial_points_list[interval_to_split_idx + 1])) break;

            Point new_point = generateNewTrialPoint(interval_to_split_idx);
            bool already_exists_nearby = false;
            for (const auto& existing_p : trial_points_list) {
                if (std::abs(existing_p.x_param - new_point.x_param) < 1e-10) {
                    already_exists_nearby = true;
                    break;
                }
            }
            if (!already_exists_nearby) {
                auto it_insert = std::lower_bound(trial_points_list.begin(), trial_points_list.end(), new_point,[](const Point& p1, const Point& p2) {
                    return p1.x_param < p2.x_param;
                });
                trial_points_list.insert(it_insert, new_point);
            }
        }

        auto min_trial_it = std::min_element(trial_points_list.begin(), trial_points_list.end(),[](const Point& t1, const Point& t2) {
            return t1.z_value < t2.z_value;
        });

        if (min_trial_it != trial_points_list.end()) return min_trial_it->y_coords;
        return {};
    }

    int GetIterationCount() const { return iteration_counter; }
    int getPointsCount() const { return static_cast<int>(trial_points_list.size()); }

    void saveGrishaginPoints(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out.is_open()) return;
        for (const auto& p : trial_points_list) {
            for (size_t i = 0; i < p.y_coords.size(); ++i) {
                out << p.y_coords[i] << (i + 1 == p.y_coords.size() ? " " : ",");
            }
            out << p.z_value << std::endl;
        }
        out.close();
    }
};

class ConstrainedMinimizer {
private:
    const IConstrainedProblem& problem;
    double a, b;
    double epsilon, r_param;
    std::vector<double> reserves_eps;
    int max_iterations;
    std::ostream& logger;
    int problem_dim, peano_m, peano_key;
    LocalizationSettings loc_settings;
    double inv_dim;

    struct TrialPoint {
        double x;
        int v;
        double z;
        std::vector<double> y_mapped;
    };
    std::vector<TrialPoint> trial_points;
    int iteration_count;

    TrialPoint performTrial(double x_val) {
        if (g_strict_budget_mode) {
            g_eval_counter++;
        }
        int m = problem.GetConstraintsCount();
        std::vector<double> y_args;
        if (problem_dim > 1) y_args = peanoMapping(x_val, peano_m, problem_dim, peano_key);
        else y_args.push_back(x_val);

        for (int j = 0; j < m; ++j) {
            double g_val = problem.ComputeConstraint(j, y_args);
            if (g_val > 0) {
                if (g_strict_budget_mode && g_debug_convergence_file.is_open()) {
                    if (g_current_best_val == std::numeric_limits<double>::infinity()) {
                        g_debug_convergence_file << "итерация " << g_eval_counter << ": N/A\n";
                    } else {
                        g_debug_convergence_file << "итерация " << g_eval_counter << ": " << std::fixed << std::setprecision(8) << g_current_best_val << "\n";
                    }
                }
                return { x_val, j + 1, g_val, y_args };
            }
        }
        double phi_val = problem.ComputeObjective(y_args);

        if (g_strict_budget_mode) {
            if (phi_val < g_current_best_val) {
                g_current_best_val = phi_val;
            }
            if (g_debug_convergence_file.is_open()) {
                g_debug_convergence_file << "итерация " << g_eval_counter << ": " << std::fixed << std::setprecision(8) << g_current_best_val << "\n";
            }
        }

        return { x_val, m + 1, phi_val, y_args };
    }

    double getDelta(double x1, double x2) const {
        double dx = std::abs(x1 - x2);
        if (problem_dim > 1) return std::pow(dx, inv_dim);
        return dx;
    }

public:
    ConstrainedMinimizer(const IConstrainedProblem& prob, double start_a, double end_b,
        double eps, double r, const std::vector<double>& reserves,
        int max_iter, std::ostream& log_stream,
        int dim = 1, int p_m = 12, int p_key = 1, LocalizationSettings loc = {})
        : problem(prob), a(start_a), b(end_b), epsilon(eps), r_param(r),
        reserves_eps(reserves), max_iterations(max_iter), logger(log_stream), iteration_count(0),
        problem_dim(dim), peano_m(p_m), peano_key(p_key), loc_settings(loc)
    {
        inv_dim = problem_dim > 0 ? 1.0 / static_cast<double>(problem_dim) : 1.0;
    }

   double findMinimum() {
       trial_points.clear();
       trial_points.reserve(max_iterations + 3);
       iteration_count = 0;

       trial_points.push_back(performTrial(a));
       trial_points.push_back(performTrial(b));
       std::sort(trial_points.begin(), trial_points.end(),[](const auto& p1, const auto& p2) { return p1.x < p2.x; });

       std::map<int, double> global_mu_estimates;
       if (trial_points[0].v == trial_points[1].v) {
           double delta = getDelta(trial_points[0].x, trial_points[1].x);
           if (delta > 1e-10) global_mu_estimates[trial_points[0].v] = std::abs(trial_points[0].z - trial_points[1].z) / delta;
       }

       std::vector<double> known_optimum_y;
       try { known_optimum_y = problem.GetOptimumPoint(); } catch (...) {}

       bool isKnownProblem = !known_optimum_y.empty() && (known_optimum_y.size() == (size_t)problem_dim);
       double effectiveEps = isKnownProblem ? (epsilon / 2.0) : epsilon;

       for (int k = 0; k < max_iterations; ++k) {
           iteration_count = k + 1;

           // Отключаем ранний выход при достижении оптимума, если задан строгий бюджет
           if (isKnownProblem && !g_strict_budget_mode) {
               double min_z_feasible = std::numeric_limits<double>::infinity();
               const TrialPoint* best_feasible = nullptr;
               int m = problem.GetConstraintsCount();
               for (const auto& p : trial_points) {
                   if (p.v == m + 1 && p.z < min_z_feasible) {
                       min_z_feasible = p.z;
                       best_feasible = &p;
                   }
               }
               if (best_feasible) {
                   bool all_coords_close = true;
                   for (int j = 0; j < problem_dim; ++j) {
                       if (std::abs(best_feasible->y_mapped[j] - known_optimum_y[j]) > effectiveEps) {
                           all_coords_close = false; break;
                       }
                   }
                   if (all_coords_close) break;
               }
           }

           std::map<int, double> mu_estimates = global_mu_estimates;
           for (int v_idx = 1; v_idx <= problem.GetConstraintsCount() + 1; ++v_idx) {
               if (mu_estimates.find(v_idx) == mu_estimates.end() || mu_estimates[v_idx] < 1e-9) mu_estimates[v_idx] = 1.0;
           }

           std::map<int, double> min_z_for_v;
           for (const auto& p : trial_points) {
               if (min_z_for_v.find(p.v) == min_z_for_v.end() || p.z < min_z_for_v[p.v]) min_z_for_v[p.v] = p.z;
           }

           std::map<int, double> z_star_estimates;
           int m = problem.GetConstraintsCount();
           for (int v = 1; v <= m + 1; ++v) {
               if (min_z_for_v.find(v) != min_z_for_v.end()) {
                   if (v <= m) z_star_estimates[v] = (min_z_for_v[v] <= 0) ? -reserves_eps[v - 1] : min_z_for_v[v];
                   else z_star_estimates[v] = min_z_for_v[v];
               }
           }

           bool is_local_step = false;
           if (loc_settings.enabled && iteration_count > loc_settings.delay) {
               int cycle = loc_settings.global_ratio + loc_settings.local_ratio;
               if (((iteration_count - loc_settings.delay - 1) % cycle) >= loc_settings.global_ratio) is_local_step = true;
           }

           double max_R = -std::numeric_limits<double>::infinity();
           size_t max_R_idx = 0;

           for (size_t i = 0; i < trial_points.size() - 1; ++i) {
               const auto& p_left = trial_points[i];
               const auto& p_right = trial_points[i + 1];
               double delta_i = getDelta(p_left.x, p_right.x);
               double R = 0;
               int v_left = p_left.v, v_right = p_right.v;

               if (v_left == v_right) {
                   int v = v_left;
                   double r_mu_v = r_param * mu_estimates[v];
                   double dz = p_right.z - p_left.z;
                   double z_star = z_star_estimates[v];
                   R = delta_i + (dz * dz) / (r_mu_v * r_mu_v * delta_i) - 2.0 * (p_left.z + p_right.z - 2.0 * z_star) / r_mu_v;
                   if (is_local_step) {
                       double diffL = p_left.z - z_star;
                       double diffR = p_right.z - z_star;
                       double delta_loc = std::sqrt(std::max(diffL, 0.0) * std::max(diffR, 0.0)) / mu_estimates[v];
                       R *= 1.0 / (delta_loc + std::pow(1.5, -loc_settings.alpha));
                   }
               } else if (v_left < v_right) {
                   R = 2 * delta_i - 4 * (p_right.z - z_star_estimates[v_right]) / (r_param * mu_estimates[v_right]);
               } else {
                   R = 2 * delta_i - 4 * (p_left.z - z_star_estimates[v_left]) / (r_param * mu_estimates[v_left]);
               }

               if (R > max_R) { max_R = R; max_R_idx = i; }
           }

           // КРИТЕРИЙ №3: Проверка длины максимального отрезка (отключена в режиме строгого бюджета)
           if (!g_strict_budget_mode) {
               double currentDelta = getDelta(trial_points[max_R_idx].x, trial_points[max_R_idx + 1].x);
               if (!isKnownProblem) { if (currentDelta <= epsilon) break; }
               else { if (currentDelta <= 1e-14) break; }
           }

           const auto& p_left = trial_points[max_R_idx];
           const auto& p_right = trial_points[max_R_idx + 1];
           double new_x;

           if (p_left.v != p_right.v) new_x = (p_left.x + p_right.x) / 2.0;
           else {
               int v = p_left.v;
               double mu_v = mu_estimates[v];
               if (problem_dim > 1) {
                   double term = std::abs(p_right.z - p_left.z) / mu_v;
                   new_x = (p_left.x + p_right.x) / 2.0 - sgn(p_right.z - p_left.z) * (1.0 / (2.0 * r_param)) * std::pow(term, (double)problem_dim);
               } else {
                   new_x = (p_left.x + p_right.x) / 2.0 - (p_right.z - p_left.z) / (2.0 * r_param * mu_v);
               }
           }

           new_x = std::max(p_left.x + 1e-12, std::min(p_right.x - 1e-12, new_x));
           TrialPoint new_point = performTrial(new_x);

           for (const auto& p : trial_points) {
               if (p.v == new_point.v) {
                   double delta = getDelta(p.x, new_point.x);
                   if (delta > 1e-10) {
                       double current_mu = std::abs(p.z - new_point.z) / delta;
                       if (global_mu_estimates.find(new_point.v) == global_mu_estimates.end() || current_mu > global_mu_estimates[new_point.v]) {
                           global_mu_estimates[new_point.v] = current_mu;
                       }
                   }
               }
           }

           auto it_insert = std::lower_bound(trial_points.begin(), trial_points.end(), new_point,[](const TrialPoint& p1, const TrialPoint& p2) {
               return p1.x < p2.x;
           });
           trial_points.insert(it_insert, new_point);
       }

       double min_phi = std::numeric_limits<double>::infinity();
       double best_x = trial_points[0].x;
       int m = problem.GetConstraintsCount();

       for (const auto& p : trial_points) {
           if (p.v == m + 1 && p.z < min_phi) {
               min_phi = p.z;
               best_x = p.x;
           }
       }
       return best_x;
   }

    int getIterationCount() const { return iteration_count; }
    int getPointsCount() const { return static_cast<int>(trial_points.size()); }

    void savePoints(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out.is_open()) return;
        for (const auto& p : trial_points) {
            for (size_t i = 0; i < p.y_mapped.size(); ++i) out << p.y_mapped[i] << (i == p.y_mapped.size() - 1 ? "" : " ");
            out << " " << p.z << " " << p.v << std::endl;
        }
        out.close();
    }
};

class GKLS_Adapter : public IConstrainedProblem {
private:
    IConstrainedOptProblem* gkls_prob;
    int dim;
    std::vector<double> lb, ub;
public:
    GKLS_Adapter(IConstrainedOptProblem* prob) : gkls_prob(prob), dim(prob->GetDimension()) {
        gkls_prob->GetBounds(lb, ub);
    }
    std::vector<double> scale(const std::vector<double>& x_unit) const {
        std::vector<double> x_real(dim);
        for (int i = 0; i < dim; ++i) x_real[i] = lb[i] + x_unit[i] * (ub[i] - lb[i]);
        return x_real;
    }
    double ComputeObjective(const std::vector<double>& x_unit) const override {
        return gkls_prob->ComputeFunction(scale(x_unit));
    }
    double ComputeConstraint(int constraint_idx, const std::vector<double>& x_unit) const override {
        return gkls_prob->ComputeConstraint(constraint_idx, scale(x_unit));
    }
    int GetConstraintsCount() const override { return gkls_prob->GetConstraintsNumber(); }
    std::vector<double> GetOptimumPoint() const override {
        std::vector<double> real_opt;
        try { real_opt = gkls_prob->GetOptimumPoint(); } catch (...) {}
        if (real_opt.empty()) return {};
        std::vector<double> unit_opt(dim);
        for (int i = 0; i < dim; ++i) unit_opt[i] = (real_opt[i] - lb[i]) / (ub[i] - lb[i]);
        return unit_opt;
    }
};

int main() {
    IOptProblem* lastGrishaginFunc = nullptr;
    double epsilon_param = 0.0, r_param = 0.0;
    int max_iterations_val = 15000;

    std::cout << "Enter r parametr:" << std::endl;
    std::cin >> r_param;
    std::cout << "Enter epsilon:" << std::endl;
    std::cin >> epsilon_param;

    std::cout << "Enter max iterations (or max evals for DIRECT, default 15000, 1000 for strict tests): ";
    std::cin >> max_iterations_val;

    std::cout << "Strict budget mode (disable early stopping & output debug convergence)? (1 - Yes, 0 - No): ";
    int strict_budget_input = 0;
    std::cin >> strict_budget_input;
    g_strict_budget_mode = (strict_budget_input == 1);

    if (g_strict_budget_mode) {
        g_debug_convergence_file.open("debug_convergence.txt");
    }

    int algorithmChoice = 0;
    std::cout << "\nChoose optimization algorithm:" << std::endl;
    std::cout << "1 - Strongin's method" << std::endl;
    std::cout << "2 - DIRECT algorithm" << std::endl;
    std::cout << "Enter your choice: ";
    std::cin >> algorithmChoice;

    LocalizationSettings loc;
    if (algorithmChoice == 1) {
        int stronginVariant = 0;
        std::cout << "\nChoose Strongin's method variant:" << std::endl;
        std::cout << "1 - Basic algorithm" << std::endl;
        std::cout << "2 - Algorithm with localization optimization (MIML)" << std::endl;
        std::cout << "Enter your choice: ";
        std::cin >> stronginVariant;

        if (stronginVariant == 2) {
            loc.enabled = true;
            std::cout << "Enter number of initial base iterations (Delay): ";
            std::cin >> loc.delay;
            std::cout << "Enter parameter alpha (refinement strength): ";
            std::cin >> loc.alpha;
            std::cout << "Enter cycle ratio - Global steps: ";
            std::cin >> loc.global_ratio;
            std::cout << "Enter cycle ratio - Local steps: ";
            std::cin >> loc.local_ratio;
        }
    }

    int taskType = 0;
    std::cout << "\nChoose task type:" << std::endl;
    std::cout << "1 - Unconstrained Hill/Shekel" << std::endl;
    std::cout << "2 - Grishagin" << std::endl;
    std::cout << "3 - Constrained 1D Example (Example3.2)" << std::endl;
    std::cout << "4 - OptSqConstr1D" << std::endl;
    std::cout << "5 - Custom Constrained Problem (Hill/Shekel combination)" << std::endl;
    std::cout << "6 - Example 5.1 (Multidimensional Constrained)" << std::endl;
    std::cout << "7 - GKLS Multidimensional Constrained" << std::endl;
    std::cout << "Enter your choice: ";
    std::cin >> taskType;

    std::ofstream log_file("minimization_log_cpp.txt", std::ios::out);
    if (!log_file.is_open()) return 1;

    int Test_func_counter = 0;
    std::vector<int> active_tasks;
    bool silent_plots = false;

    // Ввод выборки задач и исключений (Черные/Белые списки)
    if (taskType == 1 || taskType == 2 || taskType == 5 || taskType == 7) {
        std::cout << "\nEnter the number of functions in series (or 0 to read from 'custom_tasks.txt'): ";
        int raw_count;
        std::cin >> raw_count;

        if (raw_count == 0) {
            std::ifstream custom_file("custom_tasks.txt");
            if (custom_file.is_open()) {
                int val;
                while (custom_file >> val) {
                    active_tasks.push_back(val);
                }
                custom_file.close();
                std::cout << "Loaded " << active_tasks.size() << " tasks from custom_tasks.txt.\n";
            } else {
                std::cerr << "Could not open custom_tasks.txt! Create it and put indices separated by spaces.\n";
                return 1;
            }
        } else if (raw_count == 1) {
            std::cout << "Enter the specific task index (e.g. 1-100): ";
            int sidx;
            std::cin >> sidx;
            active_tasks.push_back(sidx);
        } else {
            std::cout << "Enter indices of tasks to exclude separated by space (or press Enter/input 0 to skip):\n";
            std::string exclude_line;
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::getline(std::cin, exclude_line);
            std::vector<int> exclude_list;
            if (!exclude_line.empty() && exclude_line != "0") {
                std::istringstream iss(exclude_line);
                int val;
                while (iss >> val) exclude_list.push_back(val);
            }
            for (int i = 1; i <= raw_count; ++i) {
                if (std::find(exclude_list.begin(), exclude_list.end(), i) == exclude_list.end()) {
                    active_tasks.push_back(i);
                }
            }
        }
        Test_func_counter = static_cast<int>(active_tasks.size());

        if (taskType == 5 || taskType == 7 || taskType == 2) {
            std::cout << "Do you want to save plots silently to 'plots_output' folder? (1 - Yes, 0 - No): ";
            int sp = 0;
            std::cin >> sp;
            silent_plots = (sp == 1);
            if (silent_plots) {
                std::filesystem::create_directory("plots_output");
            }
        }
    }

    THillProblemFamily hillFamily;
    TShekelProblemFamily shekelFamily;
    TGrishaginProblemFamily grishaginFamily;

    std::ofstream dataFile("graphics.txt");
    if (!dataFile.is_open()) return 1;
    dataFile << std::fixed << std::setprecision(8);

    if (taskType == 1) {
        int subTaskType = 0;
        std::cout << "Choose sub-task type: 1 - Hill, 2 - Shekel: " << std::endl;
        std::cin >> subTaskType;
        for (size_t iter_idx = 0; iter_idx < active_tasks.size(); ++iter_idx) {
            int currentIdx = active_tasks[iter_idx];
            IOptProblem* func_base = nullptr;
            if (subTaskType == 1) {
                int hillIndex = (currentIdx - 1) % hillFamily.GetFamilySize();
                func_base = hillFamily[hillIndex];
                std::cout << "\nHill Problem " << hillIndex + 1 << std::endl;
            }
            else {
                int shekelIndex = (currentIdx - 1) % shekelFamily.GetFamilySize();
                func_base = shekelFamily[shekelIndex];
                std::cout << "\nShekel Problem " << shekelIndex + 1 << std::endl;
            }

            if (!func_base) continue;

            std::vector<double> opt_pt;
            try { opt_pt = func_base->GetOptimumPoint(); } catch (...) {}

            if (!opt_pt.empty()) {
                std::cout << "Actual minimum: " << func_base->GetOptimumValue() << " at x = " << opt_pt[0] << std::endl;
            }

            vector<double> lower_bounds, upper_bounds;
            func_base->GetBounds(lower_bounds, upper_bounds);

            // Адаптация эпсилон для DIRECT
            double max_range = 1e-12;
            for (size_t i = 0; i < lower_bounds.size(); ++i) max_range = std::max(max_range, upper_bounds[i] - lower_bounds[i]);
            double direct_eps = epsilon_param / max_range;
            if (g_strict_budget_mode) direct_eps = -1.0;

            if (algorithmChoice == 1) {
                if (subTaskType == 1) {
                    auto* func = dynamic_cast<THillProblem*>(func_base);
                    if (!func) continue;
                    Minimizer<THillProblem> problem(lower_bounds, upper_bounds, epsilon_param, r_param, *func, log_file, max_iterations_val, loc);
                    auto startTime = std::chrono::high_resolution_clock::now();
                    vector<double> calculatedMinPoint = problem.findMinimum();
                    auto endTime = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration = endTime - startTime;

                    if (!calculatedMinPoint.empty()) std::cout << "Found value: " << func->ComputeFunction(calculatedMinPoint) << " at x: " << calculatedMinPoint[0] << std::endl;
                    std::cout << "Iterations: " << problem.GetIterationCount() << "\nElapsed: " << duration.count() * 1000 << " ms\n";
                    dataFile << iter_idx << " " << problem.getPointsCount() << std::endl;
                } else {
                    auto* func = dynamic_cast<TShekelProblem*>(func_base);
                    if (!func) continue;
                    Minimizer<TShekelProblem> problem(lower_bounds, upper_bounds, epsilon_param, r_param, *func, log_file, max_iterations_val, loc);
                    auto startTime = std::chrono::high_resolution_clock::now();
                    vector<double> calculatedMinPoint = problem.findMinimum();
                    auto endTime = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration = endTime - startTime;

                    if (!calculatedMinPoint.empty()) std::cout << "Found value: " << func->ComputeFunction(calculatedMinPoint) << " at x: " << calculatedMinPoint[0] << std::endl;
                    std::cout << "Iterations: " << problem.GetIterationCount() << "\nElapsed: " << duration.count() * 1000 << " ms\n";
                    dataFile << iter_idx << " " << problem.getPointsCount() << std::endl;
                }
            } else if (algorithmChoice == 2) {
                g_current_unconstrained_problem = func_base;
                g_current_constrained_problem = nullptr;
                int dim = func_base->GetDimension();
                vector<double> x_result(dim);
                double min_f_result;
                int direct_iterations = g_strict_budget_mode ? 1000000 : max_iterations_val;
                int direct_fevals = max_iterations_val;

                setup_direct_stopping(lower_bounds, upper_bounds, epsilon_param, opt_pt);
                if (Test_func_counter == 1) open_direct_log("trial_points.txt");

                auto startTime = std::chrono::high_resolution_clock::now();
                direct_optimize(objective_function_for_direct, (void*)func_base, dim, lower_bounds.data(), upper_bounds.data(),
                    x_result.data(), &min_f_result, &direct_fevals, &direct_iterations, 0.0, 0.0,
                    1e-4, 0.0, direct_eps, 0.0, &g_direct_force_stop, DIRECT_UNKNOWN_FGLOBAL, 0.0, nullptr, DIRECT_GABLONSKY);
                auto endTime = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = endTime - startTime;

                close_direct_log();
                if (g_direct_force_stop && !g_direct_stopped_x.empty()) {
                    x_result = g_direct_stopped_x; min_f_result = g_direct_stopped_f;
                }

                std::cout << "DIRECT found value: " << min_f_result << "\nIters: " << direct_iterations << "\nEvals: " << direct_fevals << "\nElapsed: " << duration.count() * 1000 << " ms\n";
                dataFile << iter_idx << " " << direct_fevals << std::endl;
            }
        }
    }
    else if (taskType == 2) {
        for (size_t iter_idx = 0; iter_idx < active_tasks.size(); ++iter_idx) {
            int currentIdx = active_tasks[iter_idx];
            int grishaginIndex = (currentIdx - 1) % grishaginFamily.GetFamilySize();

            IOptProblem* base_func = grishaginFamily[grishaginIndex];
            auto* func = dynamic_cast<TGrishaginProblem*>(base_func);
            if (!func) continue;
            lastGrishaginFunc = func;

            std::vector<double> opt_p;
            try { opt_p = func->GetOptimumPoint(); } catch (...) {}

            if (Test_func_counter == 1 && !opt_p.empty()) {
                std::cout << "\nGrishagin Problem #" << grishaginIndex + 1 << std::endl;
                std::cout << "Actual min: " << func->GetOptimumValue() << " at y = (" << opt_p[0] << ", " << opt_p[1] << ")" << std::endl;
            }

            if (algorithmChoice == 1) {
                Minimizer<TGrishaginProblem> minimizer(0.0, 1.0, epsilon_param, r_param, *func, 12, 2, 1, log_file, max_iterations_val, loc);
                auto startTime = std::chrono::high_resolution_clock::now();
                vector<double> y_res = minimizer.findMinimum();
                auto endTime = std::chrono::high_resolution_clock::now();
                chrono::duration<double> duration = endTime - startTime;

                if (Test_func_counter > 1) {
                    std::cout << "Problem " << currentIdx << " solved in " << minimizer.GetIterationCount() << " iters (" << minimizer.getPointsCount() << " points). Val: " << std::fixed << std::setprecision(6) << (!y_res.empty() ? func->ComputeFunction(y_res) : 0) << " (True: " << func->GetOptimumValue() << ").\n";
                } else {
                    if (!y_res.empty()) std::cout << "Found min: " << func->ComputeFunction(y_res) << " at y = (" << y_res[0] << ", " << y_res[1] << ")\n";
                    std::cout << "Iterations: " << minimizer.GetIterationCount() << "\nElapsed: " << duration.count() * 1000 << " ms\n";
                }

                dataFile << iter_idx << " " << minimizer.getPointsCount() << std::endl;
                if (Test_func_counter == 1 || silent_plots) minimizer.saveGrishaginPoints("grishagin_points.txt");

            } else if (algorithmChoice == 2) {
                g_current_unconstrained_problem = func;
                vector<double> lb, ub;
                func->GetBounds(lb, ub);

                // Адаптация эпсилон для DIRECT
                double max_r = 1e-12;
                for(size_t i=0; i<lb.size(); ++i) max_r = std::max(max_r, ub[i] - lb[i]);
                double direct_eps = epsilon_param / max_r;
                if (g_strict_budget_mode) direct_eps = -1.0;

                vector<double> x_result(2);
                double min_f_result;
                int direct_iterations = g_strict_budget_mode ? 1000000 : max_iterations_val;
                int direct_fevals = max_iterations_val;

                setup_direct_stopping(lb, ub, epsilon_param, opt_p);
                if (Test_func_counter == 1 || silent_plots) open_direct_log("grishagin_points.txt");

                auto startTime = std::chrono::high_resolution_clock::now();
                direct_optimize(objective_function_for_direct, nullptr, 2, lb.data(), ub.data(),
                    x_result.data(), &min_f_result, &direct_fevals, &direct_iterations, 0.0, 0.0,
                    1e-4, 0.0, direct_eps, 0.0, &g_direct_force_stop, DIRECT_UNKNOWN_FGLOBAL, 0.0, nullptr, DIRECT_GABLONSKY);
                auto endTime = std::chrono::high_resolution_clock::now();
                chrono::duration<double> duration = endTime - startTime;

                close_direct_log();
                if (g_direct_force_stop && !g_direct_stopped_x.empty()) {
                    x_result = g_direct_stopped_x; min_f_result = g_direct_stopped_f;
                }

                if (Test_func_counter > 1) {
                    std::cout << "Problem " << currentIdx << " solved in " << direct_iterations << " iters (" << direct_fevals << " points). Val: " << std::fixed << std::setprecision(6) << min_f_result << " (True: " << func->GetOptimumValue() << ").\n";
                } else {
                    std::cout << "DIRECT found min: " << min_f_result << " at y = (" << x_result[0] << ", " << x_result[1] << ")\n";
                    std::cout << "Iters: " << direct_iterations << "\nEvals: " << direct_fevals << "\nElapsed: " << duration.count() * 1000 << " ms\n";
                }
                dataFile << iter_idx << " " << direct_fevals << std::endl;
            }

            if ((Test_func_counter == 1 || silent_plots) && !opt_p.empty()) {
                string cmd = "python \"../../grishagin_visualization.py\" " + to_string(opt_p[0]) + " " + to_string(opt_p[1]);
                if (silent_plots) cmd += " silent " + to_string(currentIdx);
                system(cmd.c_str());
            }
        }
    }
    else if (taskType == 3 || taskType == 4 || taskType == 6) {
        std::string problemName;
        if (taskType == 3) problemName = "Example3.2";
        else if (taskType == 4) problemName = "OptSqConstr1D";
        else problemName = "Example5.1";

        std::cout << "\n--- Testing " << problemName << " Problem ---" << std::endl;
        unique_ptr<IConstrainedProblem> problem_ptr = (taskType == 6) ? std::make_unique<Example51Problem>() : ProblemFactory(problemName);
        if (!problem_ptr) return 1;

        double a = (taskType == 3) ? -0.6 : (taskType == 6 ? 0.0 : -1.0);
        double b = (taskType == 3) ? 2.2 : 1.0;
        int dim = (taskType == 6) ? 2 : 1;

        double found_x_coord = std::numeric_limits<double>::quiet_NaN();
        double found_min_val = 0.0;
        bool solution_found = false;
        vector<double> found_point_nd(dim);

        std::vector<double> opt_pt;
        try { opt_pt = problem_ptr->GetOptimumPoint(); } catch (...) {}

        if (algorithmChoice == 1) {
            std::vector<double> reserves(problem_ptr->GetConstraintsCount(), 0.001);
            ConstrainedMinimizer solver(*problem_ptr, a, b, epsilon_param, r_param, reserves, max_iterations_val, log_file, dim, 12, 1, loc);

            auto startTime = std::chrono::high_resolution_clock::now();
            double found_x = solver.findMinimum();
            auto endTime = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = endTime - startTime;

            solver.savePoints("trial_points.txt");

            if (!std::isnan(found_x)) {
                found_x_coord = found_x;
                if (dim > 1) {
                    found_point_nd = peanoMapping(found_x, 12, dim, 1);
                    found_min_val = problem_ptr->ComputeObjective(found_point_nd);
                } else {
                    found_point_nd = {found_x};
                    found_min_val = problem_ptr->ComputeObjective(found_point_nd);
                }
                solution_found = true;

                if (dim > 1) {
                    std::cout << std::fixed << std::setprecision(8) << "\nFound min: " << found_min_val << " at y =[" << found_point_nd[0]*4.0 << ", " << -1.0 + found_point_nd[1]*4.0 << "]\n";
                } else {
                    std::cout << std::fixed << std::setprecision(8) << "Found min: " << found_min_val << " at x = " << found_x << "\n";
                }
            } else { std::cout << "Feasible minimum not found.\n"; }
            std::cout << "Iters: " << solver.getIterationCount() << "\nElapsed: " << duration.count() * 1000 << " ms\n";
        } else if (algorithmChoice == 2) {
            std::cout << "Enter penalty coefficient for DIRECT: ";
            std::cin >> g_penalty_coefficient;
            g_current_generic_constrained_problem = problem_ptr.get();

            vector<double> lower_b(dim, a), upper_b(dim, b), x_result(dim);

            // Адаптация эпсилон для DIRECT
            double max_r = 1e-12;
            for(int d=0; d<dim; ++d) max_r = std::max(max_r, upper_b[d] - lower_b[d]);
            double direct_eps = epsilon_param / max_r;
            if (g_strict_budget_mode) direct_eps = -1.0;

            double min_f_result;
            int direct_iterations = g_strict_budget_mode ? 1000000 : max_iterations_val;
            int direct_fevals = max_iterations_val;

            setup_direct_stopping(lower_b, upper_b, epsilon_param, opt_pt);
            open_direct_log("trial_points.txt");

            auto startTime = std::chrono::high_resolution_clock::now();
            direct_optimize(generic_penalty_function_for_direct, nullptr, dim, lower_b.data(), upper_b.data(),
                x_result.data(), &min_f_result, &direct_fevals, &direct_iterations, 0.0, 0.0,
                1e-4, 0.0, direct_eps, 0.0, &g_direct_force_stop, DIRECT_UNKNOWN_FGLOBAL, 0.0, nullptr, DIRECT_GABLONSKY);
            auto endTime = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = endTime - startTime;

            close_direct_log();
            if (g_direct_force_stop && !g_direct_stopped_x.empty()) {
                x_result = g_direct_stopped_x; min_f_result = g_direct_stopped_f;
            }

            found_x_coord = x_result[0];
            found_point_nd = x_result;
            found_min_val = problem_ptr->ComputeObjective(found_point_nd);
            solution_found = true;

            std::cout << "DIRECT found min: " << min_f_result << "\nReal Obj: " << found_min_val << "\nIters: " << direct_iterations << "\nEvals: " << direct_fevals << "\nElapsed: " << duration.count() * 1000 << " ms\n";
        }

        if (solution_found) {
            if (taskType == 6) {
                std::ofstream pf2d("function_plot_data_2d.txt");
                pf2d << "3 " << found_point_nd[0]*4.0 << " " << -1.0 + found_point_nd[1]*4.0 << std::endl;
                int grid_size = 100;
                for (int i = 0; i <= grid_size; ++i) {
                    for (int j = 0; j <= grid_size; ++j) {
                        vector<double> u = {(double)i / grid_size, (double)j / grid_size};
                        pf2d << u[0] * 4.0 << " " << -1.0 + u[1] * 4.0 << " " << problem_ptr->ComputeObjective(u) << " " << problem_ptr->ComputeConstraint(0, u) << " " << problem_ptr->ComputeConstraint(1, u) << " " << problem_ptr->ComputeConstraint(2, u) << std::endl;
                    }
                }
                pf2d.close();
            }
            std::ofstream plotFile("function_plot_data.txt");
            if (plotFile.is_open()) {
                plotFile << problem_ptr->GetConstraintsCount() << " " << found_x_coord << " " << found_min_val << std::endl;
                int plot_steps = 1000;
                double step_val = (b - a) / plot_steps;
                for (int k = 0; k <= plot_steps; ++k) {
                    double x_plt = a + k * step_val;
                    std::vector<double> args = (dim > 1) ? peanoMapping(x_plt, 12, dim, 1) : vector<double>{x_plt};
                    plotFile << x_plt << " " << problem_ptr->ComputeObjective(args);
                    for (int c = 0; c < problem_ptr->GetConstraintsCount(); ++c) plotFile << " " << problem_ptr->ComputeConstraint(c, args);
                    plotFile << std::endl;
                }
                plotFile.close();
                system("python \"../../plot_1d_constrained.py\"");
                if (taskType == 6) system("python \"../../plot_2d_constrained.py\"");
            }
        }
    }
    else if (taskType == 5) {
        std::cout << "\n--- Custom Constrained Problem Series Setup ---" << std::endl;
        int objectiveType, constraintType, numConstraints;
        std::cout << "Objective type (1-Hill, 2-Shekel): "; std::cin >> objectiveType;
        std::cout << "Constraints type (1-Hill, 2-Shekel): "; std::cin >> constraintType;
        std::cout << "Number of constraints: "; std::cin >> numConstraints;

        if (algorithmChoice == 2) { std::cout << "Penalty coeff: "; std::cin >> g_penalty_coefficient; }

        for (size_t iter_idx = 0; iter_idx < active_tasks.size(); ++iter_idx) {
            int currentObjIndex = active_tasks[iter_idx];

            IOptProblem* objectiveFunc = (objectiveType == 1) ? (IOptProblem*)new THillProblem(currentObjIndex) : (IOptProblem*)new TShekelProblem(currentObjIndex);
            IConstrainedOptProblem* constrainedProblem = (constraintType == 1) ? (IConstrainedOptProblem*)new HillConstrainedProblem(objectiveFunc, numConstraints, cptNormal, 0.3, 0) : (IConstrainedOptProblem*)new TShekelConstrainedProblem(objectiveFunc, numConstraints, cptNormal, 0.3, 0);

            vector<double> lo, up;
            constrainedProblem->GetBounds(lo, up);
            int dim = constrainedProblem->GetDimension();

            std::vector<double> opt_pt;
            try { opt_pt = constrainedProblem->GetOptimumPoint(); } catch (...) {}

            int iters_taken = 0, points_used = 0;
            double found_min_val = 0, found_x_coord = 0, time_ms = 0;
            bool success = false;

            if (algorithmChoice == 1) {
                // Сброс счетчиков для дебаг файла
                g_eval_counter = 0;
                g_current_best_val = std::numeric_limits<double>::infinity();
                if (g_strict_budget_mode && g_debug_convergence_file.is_open()) {
                    g_debug_convergence_file << "\n=== Задача " << currentObjIndex << " (MIML) ===\n";
                }

                std::vector<double> reserves(numConstraints, 0.001);
                auto problem_wrapper = std::make_unique<Generic1DProblem>([constrainedProblem](double x) { return constrainedProblem->ComputeFunction({ x }); });
                for (int c = 0; c < numConstraints; ++c) problem_wrapper->AddConstraint([constrainedProblem, c](double x) { return constrainedProblem->ComputeConstraint(c, { x }); });

                ConstrainedMinimizer solver(*problem_wrapper, lo[0], up[0], epsilon_param, r_param, reserves, max_iterations_val, log_file, 1, 12, 1, loc);
                auto startTime = std::chrono::high_resolution_clock::now();
                double found_x = solver.findMinimum();
                auto endTime = std::chrono::high_resolution_clock::now();
                time_ms = std::chrono::duration<double, std::milli>(endTime - startTime).count();

                iters_taken = solver.getIterationCount(); points_used = solver.getPointsCount();
                if (!std::isnan(found_x)) { found_min_val = problem_wrapper->ComputeObjective({found_x}); found_x_coord = found_x; success = true; }
                if (Test_func_counter == 1 || silent_plots) solver.savePoints("trial_points.txt");

            } else if (algorithmChoice == 2) {
                // Сброс счетчиков для дебаг файла
                g_eval_counter = 0;
                g_current_best_val = std::numeric_limits<double>::infinity();
                if (g_strict_budget_mode && g_debug_convergence_file.is_open()) {
                    g_debug_convergence_file << "\n=== Задача " << currentObjIndex << " (DIRECT-G) ===\n";
                }

                g_current_constrained_problem = constrainedProblem;
                vector<double> x_result(dim);

                // Адаптация эпсилон для DIRECT
                double max_r = 1e-12;
                for(int d=0; d<dim; ++d) max_r = std::max(max_r, up[d] - lo[d]);
                double direct_eps = epsilon_param / max_r;
                if (g_strict_budget_mode) direct_eps = -1.0;

                double min_f_result;
                int direct_iterations = g_strict_budget_mode ? 1000000 : max_iterations_val;
                int max_evals = max_iterations_val;

                setup_direct_stopping(lo, up, epsilon_param, opt_pt);
                if (Test_func_counter == 1 || silent_plots) open_direct_log("trial_points.txt");

                auto startTime = std::chrono::high_resolution_clock::now();
                direct_optimize(penalty_function_for_direct, nullptr, dim, lo.data(), up.data(),
                                x_result.data(), &min_f_result, &max_evals, &direct_iterations,
                                0.0, 0.0, 1e-4, 0.0, direct_eps, 0.0, &g_direct_force_stop, DIRECT_UNKNOWN_FGLOBAL, 0.0, nullptr, DIRECT_GABLONSKY);
                auto endTime = std::chrono::high_resolution_clock::now();
                time_ms = std::chrono::duration<double, std::milli>(endTime - startTime).count();

                close_direct_log();
                if (g_direct_force_stop && !g_direct_stopped_x.empty()) x_result = g_direct_stopped_x;

                iters_taken = direct_iterations; points_used = max_evals;
                found_x_coord = x_result[0]; found_min_val = constrainedProblem->ComputeFunction(x_result);
                success = true;
            }

            if (Test_func_counter > 1) {
                std::cout << "Problem " << currentObjIndex << " solved in " << iters_taken << " iters (" << points_used << " pts), " << time_ms << " ms.\n";
            } else {
                std::cout << "Iterations: " << iters_taken << "\nFunc evals: " << points_used << "\nElapsed time: " << time_ms << " ms\n";
            }

            if ((Test_func_counter == 1 || silent_plots) && success) {
                std::ofstream plotFile("function_plot_data.txt");
                if (plotFile.is_open()) {
                    plotFile << constrainedProblem->GetConstraintsNumber() << " " << found_x_coord << " " << found_min_val << std::endl;
                    int plot_steps = 1000;
                    double step_val = (up[0] - lo[0]) / plot_steps;
                    for (int k = 0; k <= plot_steps; ++k) {
                        double x_plt = lo[0] + k * step_val;
                        vector<double> vec_x = {x_plt};
                        plotFile << x_plt << " " << constrainedProblem->ComputeFunction(vec_x);
                        for (int c = 0; c < constrainedProblem->GetConstraintsNumber(); ++c) plotFile << " " << constrainedProblem->ComputeConstraint(c, vec_x);
                        plotFile << std::endl;
                    }
                    plotFile.close();
                    string cmd = "python \"../../plot_1d_constrained.py\"";
                    if (silent_plots) cmd += " silent " + to_string(currentObjIndex);
                    system(cmd.c_str());
                }
            }
            dataFile << iter_idx << " " << points_used << std::endl;
            delete constrainedProblem; delete objectiveFunc;
        }
        g_current_constrained_problem = nullptr;
    }
    else if (taskType == 7) {
        std::cout << "\n    GKLS Multidimensional Constrained Setup" << std::endl;
        int gkls_class_choice, num_constrs; double gkls_fraction; int dim = 2;
        std::cout << "Choose GKLS Class (1 - Simple, 2 - Hard): "; std::cin >> gkls_class_choice;
        GKLSClass cls = (gkls_class_choice == 2) ? Hard : Simple;
        std::cout << "Enter number of constraints (e.g. 2): "; std::cin >> num_constrs;
        std::cout << "Enter feasible area fraction (0.0 to 1.0, e.g. 0.3): "; std::cin >> gkls_fraction;

        int min_pos_choice;
        std::cout << "Choose global minimum location:\n 1 - Inside feasible domain (cptInFeasibleDomain)\n 2 - On the feasible border (cptOnFeasibleBorder)\n 3 - Outside feasible domain (cptOutFeasibleDomain)\nEnter choice: ";
        std::cin >> min_pos_choice;
        EConstrainedProblemType selectedProblemType = (min_pos_choice == 2) ? cptOnFeasibleBorder : ((min_pos_choice == 3) ? cptOutFeasibleDomain : cptInFeasibleDomain);

        if (algorithmChoice == 2) { std::cout << "Enter penalty coefficient for DIRECT: "; std::cin >> g_penalty_coefficient; }

        for (size_t iter_idx = 0; iter_idx < active_tasks.size(); ++iter_idx) {
            int currentIdx = active_tasks[iter_idx];
            IConstrainedOptProblem* gkls_original = new TGKLSConstrainedProblem(selectedProblemType, gkls_fraction, 0, currentIdx, dim, cls, TD, num_constrs);

            double true_opt_val = 0.0;
            vector<double> y_real, y_true;

            try { true_opt_val = gkls_original->GetOptimumValue(); } catch (...) {}
            try { y_true = gkls_original->GetOptimumPoint(); } catch (...) {}

            int iters_taken = 0, points_used = 0; double time_ms = 0.0, found_min_val = 0.0;

            if (algorithmChoice == 1) {
                // Сброс счетчиков для дебаг файла
                g_eval_counter = 0;
                g_current_best_val = std::numeric_limits<double>::infinity();
                if (g_strict_budget_mode && g_debug_convergence_file.is_open()) {
                    g_debug_convergence_file << "\n=== Задача " << currentIdx << " (MIML) ===\n";
                }

                GKLS_Adapter adapter(gkls_original);
                std::vector<double> reserves(gkls_original->GetConstraintsNumber(), 0.001);
                ConstrainedMinimizer solver(adapter, 0.0, 1.0, epsilon_param, r_param, reserves, max_iterations_val, log_file, dim, 12, 1, loc);

                auto startTime = std::chrono::high_resolution_clock::now();
                double found_x_peano = solver.findMinimum();
                auto endTime = std::chrono::high_resolution_clock::now();
                time_ms = std::chrono::duration<double, std::milli>(endTime - startTime).count();

                vector<double> y_unit = peanoMapping(found_x_peano, 12, dim, 1);
                found_min_val = adapter.ComputeObjective(y_unit);
                y_real = adapter.scale(y_unit);
                iters_taken = solver.getIterationCount(); points_used = solver.getPointsCount();

                if (Test_func_counter == 1 || silent_plots) solver.savePoints("trial_points.txt");

            } else if (algorithmChoice == 2) {
                // Сброс счетчиков для дебаг файла
                g_eval_counter = 0;
                g_current_best_val = std::numeric_limits<double>::infinity();
                if (g_strict_budget_mode && g_debug_convergence_file.is_open()) {
                    g_debug_convergence_file << "\n=== Задача " << currentIdx << " (DIRECT-G) ===\n";
                }

                g_current_constrained_problem = gkls_original;
                vector<double> lb, ub; gkls_original->GetBounds(lb, ub); y_real.resize(dim);

                // Адаптация эпсилон для DIRECT
                double max_r = 1e-12;
                for(int d=0; d<dim; ++d) max_r = std::max(max_r, ub[d] - lb[d]);
                double direct_eps = epsilon_param / max_r;
                if (g_strict_budget_mode) direct_eps = -1.0;

                double min_f_result;
                int direct_iterations = g_strict_budget_mode ? 1000000 : max_iterations_val;
                int max_evals = max_iterations_val;

                setup_direct_stopping(lb, ub, epsilon_param, y_true);
                if (Test_func_counter == 1 || silent_plots) open_direct_log("trial_points.txt");

                auto startTime = std::chrono::high_resolution_clock::now();
                direct_optimize(penalty_function_for_direct, nullptr, dim, lb.data(), ub.data(),
                    y_real.data(), &min_f_result, &max_evals, &direct_iterations, 0.0, 0.0, 1e-4, 0.0, direct_eps, 0.0, &g_direct_force_stop, DIRECT_UNKNOWN_FGLOBAL, 0.0, nullptr, DIRECT_GABLONSKY);
                auto endTime = std::chrono::high_resolution_clock::now();
                time_ms = std::chrono::duration<double, std::milli>(endTime - startTime).count();

                close_direct_log();
                if (g_direct_force_stop && !g_direct_stopped_x.empty()) { y_real = g_direct_stopped_x; }
                found_min_val = gkls_original->ComputeFunction(y_real);
                points_used = max_evals;
                iters_taken = direct_iterations;
            }

            if (Test_func_counter > 1) {
                std::cout << "Problem " << currentIdx << " solved in " << iters_taken << " iters (" << points_used << " points). Val: " << std::fixed << std::setprecision(6) << found_min_val << " (True: " << true_opt_val << ").\n";
            } else {
                std::cout << "\n    Solving GKLS Problem #" << currentIdx << "\n";
                std::cout << "Found value:  " << std::fixed << std::setprecision(8) << found_min_val << "\nTrue value:   " << true_opt_val << "\n";
                std::cout << "Found point:[" << y_real[0] << ", " << y_real[1] << "]\n";
                if (!y_true.empty()) std::cout << "True point:[" << y_true[0] << ", " << y_true[1] << "]\n";
                std::cout << "Iters: " << iters_taken << "\nEvals: " << points_used << "\nElapsed: " << time_ms << " ms\n";
            }

            if (dim == 2 && (Test_func_counter == 1 || silent_plots) && !y_true.empty()) {
                std::ofstream pf2d("function_plot_data_2d.txt");
                if (pf2d.is_open()) {
                    pf2d << gkls_original->GetConstraintsNumber() << " " << y_true[0] << " " << y_true[1] << std::endl;
                    int grid_size = 100;
                    for (int ix = 0; ix <= grid_size; ++ix) {
                        for (int iy = 0; iy <= grid_size; iy++) {
                            vector<double> u = { (double)ix / grid_size, (double)iy / grid_size };
                            GKLS_Adapter adapter(gkls_original); vector<double> r = adapter.scale(u);
                            pf2d << r[0] << " " << r[1] << " " << adapter.ComputeObjective(u);
                            for (int c = 0; c < gkls_original->GetConstraintsNumber(); ++c) pf2d << " " << adapter.ComputeConstraint(c, u);
                            pf2d << std::endl;
                        }
                    }
                    pf2d.close();
                    string cmd = "python \"../../plot_gkls_2d.py\"";
                    if (silent_plots) cmd += " silent " + to_string(currentIdx);
                    system(cmd.c_str());
                }
            }
            dataFile << iter_idx << " " << points_used << std::endl;
            delete gkls_original;
        }
    }

    std::ofstream statsFile("stats.txt");
    statsFile << epsilon_param << " " << r_param << " " << Test_func_counter << std::endl;
    statsFile.close(); dataFile.close(); log_file.close();

    if (g_strict_budget_mode) {
        g_debug_convergence_file.close();
        std::cout << "\nConvergence log successfully written to debug_convergence.txt\n";
    }

    if ((taskType == 1 || taskType == 2 || taskType == 5 || taskType == 7) && Test_func_counter > 1) system("python \"../../graphics.py\"");

    return 0;
}
