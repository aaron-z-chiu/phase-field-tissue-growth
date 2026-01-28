#ifndef PHASE_FIELD_MODEL_H
#define PHASE_FIELD_MODEL_H

#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>


class PhaseFieldModel {
public:
    struct Parameters {
        double epsilon;
        double lambda0;
        double alpha;
        double beta0;
        double dt;
        double dx;
        double dy;
        int Nx, Ny;
        int total_time_steps;
        double Lx_min, Lx_max, Ly_min, Ly_max;
        double S = 2.0; 
        double tol = 1.0e-12;
        
        Parameters() : epsilon(0.01), lambda0(10.0), alpha(0.1), beta0(0.1),
                      dt(1e-5), dx(0.01), dy(0.01), Nx(100), Ny(100),
                      total_time_steps(10000), Lx_min(0.0), Lx_max(1.0),
                      Ly_min(0.0), Ly_max(1.0) {}
    };
    

    
private:
    Parameters params;
    std::vector<std::vector<double>> phi;
    std::vector<std::vector<double>> phi_star;
    std::vector<std::vector<double>> phi_new;

    
    double M_mobility(double phi_val) const;
    double F_potential(double phi_val) const;
    double F_prime(double phi_val) const;

    double compute_same_grid_error(const std::vector<std::vector<double>> &phi1,
                                   const std::vector<std::vector<double>> &phi2) const;
    

    void solve_D1(double dt);
    void solve_D2(double dt);
    void solve_D3(double dt);

    void solve_3step(double dt_step);
    
public:
    PhaseFieldModel();
    PhaseFieldModel(const Parameters& p);
    ~PhaseFieldModel();
    
    void initialize_phi(const std::string& initial_condition_type, 
                       const std::vector<double>& condition_params = {});
    
    void solve();
    
    void run_convergence_test();
    void run_parameter_effect_test_lambda0();
    void run_parameter_effect_test_alpha();

    
    double compute_convergence_error(const std::vector<std::vector<double>>& phi_h,
                                   const std::vector<std::vector<double>>& phi_h2) const;
    
    void save_data_for_visualization(const std::string& filename, int time_step = -1) const;
    void save_phi_field(const std::string& filename) const;
    
    void set_parameters(const Parameters& p) { params = p; resize_arrays(); }
    void update_parameters_safe(const Parameters& p) { 
        Parameters old_params = params;
        params = p; 
        try {
            validate_parameters();
        } catch (const std::invalid_argument& e) {
            params = old_params;  
            throw;
        }
    } 
    const Parameters& get_parameters() const { return params; }
    const std::vector<std::vector<double>>& get_phi() const { return phi; }

    
private:
    void resize_arrays();

    void initialize_params_condition();
    void initialize_convergence_condition();
    
    void run_single_parameter_test(double lambda0_val, const std::string& output_prefix);
    void run_single_alpha_test(double alpha_val, const std::string& output_prefix);
    
    
    
    void validate_parameters() const;
};

#endif