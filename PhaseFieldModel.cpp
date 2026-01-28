#include "PhaseFieldModel.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <filesystem>
#include <sstream>

PhaseFieldModel::PhaseFieldModel()
{
    resize_arrays();
}

PhaseFieldModel::PhaseFieldModel(const Parameters &p) : params(p)
{
    resize_arrays();
    validate_parameters();
}

PhaseFieldModel::~PhaseFieldModel()
{
}

void PhaseFieldModel::resize_arrays()
{
    phi.assign(params.Nx, std::vector<double>(params.Ny, 0.0));
    phi_star.assign(params.Nx, std::vector<double>(params.Ny, 0.0));
    phi_new.assign(params.Nx, std::vector<double>(params.Ny, 0.0));
}

void PhaseFieldModel::validate_parameters() const
{
    if (params.epsilon <= 0)
    {
        throw std::invalid_argument("Invalid parameter: epsilon must be positive, got " + std::to_string(params.epsilon));
    }
    if (params.lambda0 < 0)
    {
        throw std::invalid_argument("Invalid parameter: lambda0 must be non-negative, got " + std::to_string(params.lambda0));
    }
    if (params.alpha < 0 || params.alpha > 5.0)
    {
        throw std::invalid_argument("Invalid parameter: alpha must be in [0, 5.0], got " + std::to_string(params.alpha));
    }
    if (params.beta0 < 0)
    {
        throw std::invalid_argument("Invalid parameter: beta0 must be non-negative, got " + std::to_string(params.beta0));
    }
    if (params.dt <= 0)
    {
        throw std::invalid_argument("Invalid parameter: dt must be positive, got " + std::to_string(params.dt));
    }
    if (params.dx <= 0)
    {
        throw std::invalid_argument("Invalid parameter: dx must be positive, got " + std::to_string(params.dx));
    }
    if (params.dy <= 0)
    {
        throw std::invalid_argument("Invalid parameter: dy must be positive, got " + std::to_string(params.dy));
    }
    if (params.Nx <= 0)
    {
        throw std::invalid_argument("Invalid parameter: Nx must be positive, got " + std::to_string(params.Nx));
    }
    if (params.Ny <= 0)
    {
        throw std::invalid_argument("Invalid parameter: Ny must be positive, got " + std::to_string(params.Ny));
    }
    if (params.total_time_steps <= 0)
    {
        throw std::invalid_argument("Invalid parameter: total_time_steps must be positive, got " + std::to_string(params.total_time_steps));
    }
    if (params.S < 0.0) {
        throw std::invalid_argument("Invalid parameter: S must be >= 0");
    }


    
}

double PhaseFieldModel::M_mobility(double phi_val) const
{
    return 4.0 * (phi_val - 0.5) * (phi_val - 0.5); 
}

double PhaseFieldModel::F_potential(double phi_val) const
{
    return 0.25 * phi_val * phi_val * (phi_val - 1.0) * (phi_val - 1.0);
}

double PhaseFieldModel::F_prime(double phi_val) const
{
    return phi_val * (phi_val - 1.0) * (phi_val - 0.5);
}


void PhaseFieldModel::solve_D1(double dt)
{
    const int Nx = params.Nx, Ny = params.Ny;
    const double dx = params.dx, eps2 = params.epsilon * params.epsilon;
    const double S  = params.S;

    
    const double tol = params.tol; 
    const int max_inner_iter = 1000; 

    
    phi_star = phi; 

    
    for (int iter = 0; iter < max_inner_iter; ++iter) {
        double max_diff = 0.0;

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                const double Mij   = M_mobility(phi[i][j]); 
                const double Fpij  = F_prime(phi[i][j]);
                const double ax    = (dt * Mij) / (dx * dx);
                const double sigma = dt * Mij * (S / eps2);

                
                const double RHS = phi[i][j]
                                 + dt * Mij * ( -Fpij / (2.0 * eps2) )
                                 + sigma * phi[i][j];

                
                
                double left  = (i > 0)      ? phi_star[i-1][j] : phi_star[i][j]; 
                double right = (i+1 < Nx)   ? phi_star[i+1][j] : phi_star[i][j];

                
                
                const double denom = 1.0 + 2.0 * ax + sigma;
                const double new_val = (RHS + ax * (left + right)) / denom;

                
                double diff = std::abs(new_val - phi_star[i][j]);
                if (diff > max_diff) max_diff = diff;

                
                phi_star[i][j] = new_val;
            }
        }

        
        if (max_diff < tol) {
            
            break;
        }
    }
}


void PhaseFieldModel::solve_D2(double dt)
{
    const int Nx = params.Nx, Ny = params.Ny;
    const double dy = params.dy, eps2 = params.epsilon * params.epsilon;
    const double S  = params.S;

    const double tol = params.tol;
    const int max_inner_iter = 1000;

    
    phi_new = phi_star;

    for (int iter = 0; iter < max_inner_iter; ++iter) {
        double max_diff = 0.0;

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                const double Mij   = M_mobility(phi_star[i][j]); 
                const double Fpij  = F_prime(phi_star[i][j]);
                const double ay    = (dt * Mij) / (dy * dy);
                const double sigma = dt * Mij * (S / eps2);

                const double RHS = phi_star[i][j]
                                 + dt * Mij * ( -Fpij / (2.0 * eps2) )
                                 + sigma * phi_star[i][j];

                double down = (j > 0)     ? phi_new[i][j-1] : phi_new[i][j];
                double up   = (j+1 < Ny)  ? phi_new[i][j+1] : phi_new[i][j];

                const double denom = 1.0 + 2.0 * ay + sigma;
                const double new_val = (RHS + ay * (down + up)) / denom;

                double diff = std::abs(new_val - phi_new[i][j]);
                if (diff > max_diff) max_diff = diff;

                phi_new[i][j] = new_val;
            }
        }

        if (max_diff < tol) {
            break;
        }
    }
}


void PhaseFieldModel::solve_D3(double dt)
{
    
    
    const int Nx = params.Nx, Ny = params.Ny;
    const double dx = params.dx, dy = params.dy;
    const double beta0 = params.beta0, lam0 = params.lambda0, alpha = params.alpha;
    const double eps_grad = 1e-12;

    auto clamp_idx = [&](int ii, int lo, int hi) {
        return (ii < lo) ? lo : (ii > hi) ? hi : ii;
    };
    auto val_at = [&](int ii, int jj)->double {
        ii = clamp_idx(ii, 0, Nx-1);
        jj = clamp_idx(jj, 0, Ny-1);
        return phi_new[ii][jj]; 
    };
    auto normal_x_at = [&](int ii, int jj)->double {
        const double gx_ = (val_at(ii+1,jj) - val_at(ii-1,jj)) / (2.0 * dx);
        const double gy_ = (val_at(ii,jj+1) - val_at(ii,jj-1)) / (2.0 * dy);
        const double g_  = std::sqrt(gx_*gx_ + gy_*gy_ + eps_grad);
        return gx_ / g_;
    };
    auto normal_y_at = [&](int ii, int jj)->double {
        const double gx_ = (val_at(ii+1,jj) - val_at(ii-1,jj)) / (2.0 * dx);
        const double gy_ = (val_at(ii,jj+1) - val_at(ii,jj-1)) / (2.0 * dy);
        const double g_  = std::sqrt(gx_*gx_ + gy_*gy_ + eps_grad);
        return gy_ / g_;
    };

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            
            const double nx_r = normal_x_at(i+1,j);
            const double nx_l = normal_x_at(i-1,j);
            const double ny_u = normal_y_at(i,j+1);
            const double ny_d = normal_y_at(i,j-1);
            const double dnx_dx = (nx_r - nx_l) / (2.0 * dx);
            const double dny_dy = (ny_u - ny_d) / (2.0 * dy);
            const double kappa  = dnx_dx + dny_dy;

            
            double hat = phi_new[i][j];
            if (hat < 1e-8) hat = 1e-8; else if (hat > 1.0 - 1e-8) hat = 1.0 - 1e-8;

            const double expo = -lam0 * dt * (1.0 + alpha * std::tanh(beta0 * hat * (1.0 - hat) * kappa) );
            const double g    = std::exp(expo);
            const double num  = hat;
            const double den  = hat + (1.0 - hat) * g;

            phi[i][j] = num / den;   
        }
    }
}

void PhaseFieldModel::solve_3step(double dt_step)
{
    
    solve_D1(dt_step); 
    solve_D2(dt_step); 
    solve_D3(dt_step); 
}

void PhaseFieldModel::solve()
{
    auto t0 = std::chrono::high_resolution_clock::now();

    const int max_it = params.total_time_steps;
    std::cout << "Max iterations: " << max_it << std::endl;

    std::vector<std::vector<double>> phi_n(params.Nx, std::vector<double>(params.Ny));
    std::vector<std::vector<double>> phi_dt1(params.Nx, std::vector<double>(params.Ny));

    for (int it = 1; it <= max_it; ++it) {

        
        phi_n = phi;
        solve_3step(params.dt);
        phi_dt1 = phi;

        
        phi = phi_n;
        solve_3step(params.dt / 2.0);
        solve_3step(params.dt / 2.0);

        for (int i = 0; i < params.Nx; ++i) {
            for (int j = 0; j < params.Ny; ++j) {
                phi[i][j] = 2.0 * phi[i][j] - phi_dt1[i][j];
            }
        }

        
        if (it % 100 == 0) {
            char fname[64];
            std::snprintf(fname, sizeof(fname), "iter_%05d.dat", it);
            save_phi_field(fname);
            std::cout << "Saved " << fname << std::endl;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "Solve completed, cost " << ms << " ms" << std::endl;
}



void PhaseFieldModel::initialize_phi(const std::string &initial_condition_type,
                                     const std::vector<double> &condition_params)
{
    (void)condition_params;
    if (initial_condition_type == "convergence_test")
    {
        initialize_convergence_condition();
    }
    else if (initial_condition_type == "params_effect")
    {
        initialize_params_condition();
    }
    else
    {
        std::cerr << "No proper initialization type" << std::endl;
        return;
    }
}

void PhaseFieldModel::initialize_convergence_condition()
{
    
    for (int i = 0; i < params.Nx; ++i) {

        for (int j = 0; j < params.Ny; ++j) {

            double x = params.Lx_min + (i + 0.5) * params.dx; 
            double y = params.Ly_min + (j + 0.5) * params.dy; 

            const double numer = 0.5 - y + 0.1 * std::cos(2.0 * M_PI * x);
            const double denom = 2.0 * std::sqrt(2.0) * params.epsilon;
            phi[i][j] = 0.5 * (1.0 + std::tanh(numer / denom));
        }
    }
}

void PhaseFieldModel::initialize_params_condition()
{
    for (int i = 0; i < params.Nx; ++i)
    {
        for (int j = 0; j < params.Ny; ++j)
        {
            double x = params.Lx_min + (i + 0.5) * params.dx; 
            double y = params.Ly_min + (j + 0.5) * params.dy; 

            double numerator;

            
            if (x < 3.0)
            {
                numerator = 0.4 - y + 0.3 * std::cos(2.0 / 3.0 * M_PI * x); 
            }
            else
            {
                numerator = 0.4 - y + 0.3 * std::cos(2.0 * M_PI * x);  
            }

            
            double denominator = 2.0 * sqrt(2.0) * params.epsilon;
            double argument = numerator / denominator;

            
            phi[i][j] = 0.5 + 0.5 * tanh(argument);
        }
    }
}




double PhaseFieldModel::compute_same_grid_error(const std::vector<std::vector<double>> &phi1,
                                                const std::vector<std::vector<double>> &phi2) const
{
    if (phi1.empty() || phi2.empty()) return 0.0;
    
    int Nx = phi1.size();
    int Ny = phi1[0].size();
    
    double sum = 0.0;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double diff = phi1[i][j] - phi2[i][j];
            sum += diff * diff;
        }
    }
    return std::sqrt(sum / (Nx * Ny));
}

double PhaseFieldModel::compute_convergence_error(const std::vector<std::vector<double>> &phi_h,
    const std::vector<std::vector<double>> &phi_h2) const
{
    const int Nx_c = static_cast<int>(phi_h.size());
    const int Ny_c = Nx_c ? static_cast<int>(phi_h[0].size()) : 0;
    const int Nx_f = static_cast<int>(phi_h2.size());
    const int Ny_f = Nx_f ? static_cast<int>(phi_h2[0].size()) : 0;

    if (Nx_c == 0 || Ny_c == 0 || Nx_f == 0 || Ny_f == 0) return 0.0;

    
    
    auto clamp = [](int v, int lo, int hi){ return v < lo ? lo : (v > hi ? hi : v); };

    double accum = 0.0;
    for (int i = 0; i < Nx_c; ++i) {
        int i0 = clamp(2*i,     0, Nx_f - 1);
        int i1 = clamp(2*i + 1, 0, Nx_f - 1);
        for (int j = 0; j < Ny_c; ++j) {
            int j0 = clamp(2*j,     0, Ny_f - 1);
            int j1 = clamp(2*j + 1, 0, Ny_f - 1);

            const double avg_fine =
                0.25 * (phi_h2[i0][j0] + phi_h2[i1][j0]
                      + phi_h2[i0][j1] + phi_h2[i1][j1]);

            const double diff = phi_h[i][j] - avg_fine;
            accum += diff * diff;
        }
    }

    const double mean_sq = accum / static_cast<double>(Nx_c * Ny_c);
    return std::sqrt(mean_sq);
}



void PhaseFieldModel::run_convergence_test()
{
    std::cout << "\n=== Temporal Convergence Test (Fixed h, Varying dt) ===\n";

    
    
    const double h_fixed = std::pow(2.0, -7); 
    
    
    
    const double eps = 8.0 * h_fixed / (4.0 * std::sqrt(2.0) * std::atanh(0.9));

    
    
    
    
     double dt_stability_ref = pow(2.0,-21);
    
    
    
    double dt_start = dt_stability_ref; 
    
    
    std::vector<double> dts = { dt_start, dt_start/2.0, dt_start/4.0, dt_start/8.0, dt_start/16.0 };
    
    
    
    const double T_physical = 100.0 * dt_start; 

    std::vector<std::vector<std::vector<double>>> solutions;
    solutions.reserve(dts.size());

    std::cout << "Fixed Grid h = " << h_fixed << ", Epsilon = " << eps << "\n";
    std::cout << "Physical Time T = " << T_physical << "\n";
    std::cout << "---------------------------------------------------------\n";

    
    for (size_t k = 0; k < dts.size(); ++k) {
        Parameters p;
        p.Lx_min = 0.0; p.Lx_max = 1.0;
        p.Ly_min = 0.0; p.Ly_max = 1.0;
        
        
        p.dx = h_fixed;
        p.dy = h_fixed;
        p.Nx = static_cast<int>((p.Lx_max - p.Lx_min) / p.dx);
        p.Ny = static_cast<int>((p.Ly_max - p.Ly_min) / p.dy);

        p.epsilon = eps;
        p.alpha   = 0.1;
        p.beta0   = 0.1;
        p.lambda0 = 10.0;
        
        
        
        p.tol = 1e-12; 

        
        p.dt = dts[k];
        
        
        p.total_time_steps = static_cast<int>(std::round(T_physical / p.dt));

        PhaseFieldModel m(p);
        m.initialize_phi("convergence_test"); 

        std::cout << "Running dt[" << k << "] = " << std::scientific << p.dt 
                  << " (" << p.total_time_steps << " steps)... " << std::flush;
        
        m.solve();
        
        solutions.push_back(m.get_phi());
        std::cout << "Done." << std::endl;
    }

    std::string out_file_name = "temporal_convergence_results.txt";
    std::ofstream file(out_file_name);
    
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open " << out_file_name << " for writing." << std::endl;
    }

    
    std::cout << "\n=== Temporal Accuracy Results ===\n";
    std::cout << "dt_coarse     dt_fine       L2_Error      Rate\n";
    std::cout << "---------------------------------------------------------\n";

    if (file.is_open()) {
        file << "=== Temporal Accuracy Results ===\n";
        file << "dt_coarse     dt_fine       L2_Error      Rate\n";
        file << "---------------------------------------------------------\n";
    }

    std::vector<double> errors;
    for (size_t k = 0; k < dts.size() - 1; ++k) {
        
        
        double diff = compute_same_grid_error(solutions[k], solutions[k+1]);
        errors.push_back(diff);
    }

    for (size_t k = 0; k < errors.size(); ++k) {
        double dt_c = dts[k];
        double dt_f = dts[k+1];
        
        std::cout << std::scientific << std::setprecision(4) 
                  << dt_c << "  " << dt_f << "  " << errors[k];

        if (file.is_open()) {
            file << std::scientific << std::setprecision(4) 
                 << dt_c << "  " << dt_f << "  " << errors[k];
        }
        
        if (k > 0) {
            
            double rate = std::log2(errors[k-1] / errors[k]);
            std::cout << "  " << std::fixed << std::setprecision(2) << rate;
            if (file.is_open()) {
                file << "  " << std::fixed << std::setprecision(2) << rate;
            }
        } else {
            std::cout << "  -";
            if (file.is_open()) file << "  -";
        }
        std::cout << "\n";
    }
    std::cout << "---------------------------------------------------------\n";
    if (file.is_open()) {
        file << "---------------------------------------------------------\n";
        file.close();
        std::cout << "Results saved to " << out_file_name << std::endl;
    }
}

void PhaseFieldModel::run_single_parameter_test(double lambda0_val, const std::string &output_prefix)
{
    Parameters temp_params = params;
    temp_params.lambda0 = lambda0_val;


    temp_params.dt = 4 * temp_params.epsilon * temp_params.epsilon * temp_params.dx * temp_params.dx /
                     (temp_params.dx * temp_params.dx + 16.0 * temp_params.epsilon * temp_params.epsilon);


    std::cout << "\nTesting λ₀=" << lambda0_val << " Effect" << std::endl;
    std::cout << "Nx: " << temp_params.Nx << ", Ny" << temp_params.Ny
              << ", Timesteps: " << temp_params.total_time_steps << std::endl;

    PhaseFieldModel model(temp_params);
    model.initialize_phi("params_effect");

    
    std::string init_filename = output_prefix + "_lambda0_" + std::to_string(static_cast<int>(lambda0_val)) + "_t0.dat";
    model.save_phi_field(init_filename);

    model.solve();

    std::string filename = output_prefix + "_lambda0_" + std::to_string(static_cast<int>(lambda0_val)) + ".dat";
    model.save_phi_field(filename);

    


    std::cout << "λ₀=" << lambda0_val << " Results: " << std::endl;

}

void PhaseFieldModel::run_parameter_effect_test_lambda0()
{
    run_single_parameter_test(30.0, "fig");
    run_single_parameter_test(100.0, "fig");
}

void PhaseFieldModel::run_single_alpha_test(double alpha_val, const std::string &output_prefix)
{
    Parameters temp_params = params;
    temp_params.alpha = alpha_val;

    temp_params.dt = 4 * temp_params.epsilon * temp_params.epsilon * temp_params.dx * temp_params.dx /
                     (temp_params.dx * temp_params.dx + 16.0 * temp_params.epsilon * temp_params.epsilon);

    std::cout << "\nTesting α =" << alpha_val << " effect" << std::endl;
    std::cout << "Nx " << temp_params.Nx << "Ny" << temp_params.Ny
              << ", Time steps: " << temp_params.total_time_steps << std::endl;

    PhaseFieldModel model(temp_params);
    model.initialize_phi("params_effect");

    
    std::string init_filename = output_prefix + "_alpha_" + std::to_string(static_cast<double>(alpha_val)) + "_t0.dat";
    model.save_phi_field(init_filename);

    model.solve();

    std::string filename = output_prefix + "_alpha_" + std::to_string(static_cast<double>(alpha_val)) + ".dat";
    model.save_phi_field(filename);



  

    std::cout << "α=" << alpha_val << " results:" << std::endl;

}

void PhaseFieldModel::run_parameter_effect_test_alpha()
{
    run_single_alpha_test(0.05, "fig");
    run_single_alpha_test(0.5, "fig");
    run_single_alpha_test(1, "fig");
}


void PhaseFieldModel::save_phi_field(const std::string &filename) const
{
    const std::filesystem::path data_dir = std::filesystem::path("./data");
    std::error_code ec;
    std::filesystem::create_directories(data_dir, ec);
    const std::filesystem::path dest_path = data_dir / std::filesystem::path(filename).filename();
    std::ofstream file(dest_path.string());
    if (!file.is_open())
    {
        std::cerr << "Can not open the file: " << dest_path.string() << std::endl;
        return;
    }

    file << std::fixed << std::setprecision(8);
    file << "# Phase Field Model Output Data\n";
    file << "# Parameters:\n";
    file << "# epsilon=" << params.epsilon << " lambda0=" << params.lambda0
         << " alpha=" << params.alpha << " beta0=" << params.beta0 << "\n";
    file << "# Grid: Nx=" << params.Nx << " Ny=" << params.Ny
         << " dx=" << params.dx << " dy=" << params.dy << "\n";
    file << "# Domain: [" << params.Lx_min << "," << params.Lx_max << "] x ["
         << params.Ly_min << "," << params.Ly_max << "]\n";
    file << "# Data format: x y phi curvature growth_rate\n";

    
    for (int i = 0; i < params.Nx; ++i)
    {
        for (int j = 0; j < params.Ny; ++j)
        {
            double x = params.Lx_min + i * params.dx;
            double y = params.Ly_min + j * params.dy;
            double phi_val = phi[i][j];

            file << x << " " << y << " " << phi_val << "\n";
        }
    }

    
    double min_phi = 1e10, max_phi = -1e10, avg_phi = 0.0;
    int interface_points = 0;

    for (int i = 0; i < params.Nx; ++i)
    {
        for (int j = 0; j < params.Ny; ++j)
        {
            double phi_val = phi[i][j];
            min_phi = std::min(min_phi, phi_val);
            max_phi = std::max(max_phi, phi_val);
            avg_phi += phi_val;

        }
    }
    avg_phi /= (params.Nx * params.Ny);


    file << "# Statistics:\n";
    file << "# phi: min=" << min_phi << " max=" << max_phi << " avg=" << avg_phi << "\n";
    file << "# interface_points=" << interface_points << "\n";

    file.close();
    std::cout << "Data saved to: " << dest_path.string() << " (range of φ: [" << min_phi << ", " << max_phi
              << "], average: " << avg_phi << ")" << std::endl;
}

void PhaseFieldModel::save_data_for_visualization(const std::string &filename, int time_step) const
{
    const std::filesystem::path data_dir = std::filesystem::path("./data");
    std::error_code ec;
    std::filesystem::create_directories(data_dir, ec);
    const std::filesystem::path dest_path = data_dir / std::filesystem::path(filename).filename();
    std::ofstream file(dest_path.string());
    if (!file.is_open())
    {
        std::cerr << "Can not open the visualization file: " << dest_path.string() << std::endl;
        return;
    }

    file << std::fixed << std::setprecision(8);
    file << "# Phase field visualization data\n";
    file << "# Time step: " << time_step << "\n";
    file << "# Parameters: epsilon=" << params.epsilon << " lambda0=" << params.lambda0
         << " alpha=" << params.alpha << " beta0=" << params.beta0 << "\n";
    file << "# Domain: [" << params.Lx_min << "," << params.Lx_max << "] x ["
         << params.Ly_min << "," << params.Ly_max << "]\n";
    file << "# Grid: " << params.Nx << " x " << params.Ny << "\n";
    file << "# Format: x y phi gradient_x gradient_y gradient_magnitude curvature growth_rate mobility\n";

    
    for (int i = 0; i < params.Nx; ++i)
    {
        for (int j = 0; j < params.Ny; ++j)
        {
            double x = params.Lx_min + i * params.dx;
            double y = params.Ly_min + j * params.dy;
            double phi_val = phi[i][j];

        


            
            double mobility = M_mobility(phi_val);

            file << x << " " << y << " " << phi_val << " "
                 << mobility << "\n";
        }
        file << "\n"; 
    }

    
    file << "\n# Interface contour (phi=0.5)\n";
    file << "# x y\n";
    for (int i = 0; i < params.Nx - 1; ++i)
    {
        for (int j = 0; j < params.Ny - 1; ++j)
        {
            double phi00 = phi[i][j];
            double phi10 = phi[i + 1][j];
            double phi01 = phi[i][j + 1];

            
            double iso_value = 0.5;

            
            if ((phi00 - iso_value) * (phi10 - iso_value) < 0)
            {
                double t = (iso_value - phi00) / (phi10 - phi00);
                double x_inter = params.Lx_min + (i + t) * params.dx;
                double y_inter = params.Ly_min + j * params.dy;
                file << x_inter << " " << y_inter << "\n";
            }
            if ((phi00 - iso_value) * (phi01 - iso_value) < 0)
            {
                double t = (iso_value - phi00) / (phi01 - phi00);
                double x_inter = params.Lx_min + i * params.dx;
                double y_inter = params.Ly_min + (j + t) * params.dy;
                file << x_inter << " " << y_inter << "\n";
            }
        }
    }

    
    file << "\n# Statistics Summary\n";
    double min_phi = 1e10, max_phi = -1e10;
    double interface_area = 0.0;


    for (int i = 0; i < params.Nx; ++i)
    {
        for (int j = 0; j < params.Ny; ++j)
        {
            double phi_val = phi[i][j];
            min_phi = std::min(min_phi, phi_val);
            max_phi = std::max(max_phi, phi_val);

            
        }
    }

    file << "# phi_range: [" << min_phi << ", " << max_phi << "]\n";
    file << "# interface_area: " << interface_area << "\n";
    file << "# time_step: " << time_step << "\n";
    file << "# physical_time: " << time_step * params.dt << "\n";

    file.close();
    std::cout << "Virsualization Data Saved: " << dest_path.string()
              << " (including " << params.Nx * params.Ny << " points)" << std::endl;
}