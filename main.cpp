#include "PhaseFieldModel.h"

void run_convergence_analysis() {
    std::cout << "\n=== Start Convergence Test ===" << std::endl;    
    PhaseFieldModel m;
    m.run_convergence_test();
}

void run_lambda0_parameter_study() {
    std::cout << "\n=== λ₀ Effect ===" << std::endl;
    PhaseFieldModel::Parameters params;
    params.dx = 0.01; params.dy = 0.01;
    params.epsilon = 5 * params.dx / (4.0 * sqrt(2.0) * atanh(0.9));
    params.alpha = 0.5;
    params.beta0 = 0.1;
    params.total_time_steps = 20000;
    params.Lx_min = 0.0; params.Lx_max = 4.0;
    params.Ly_min = 0.0; params.Ly_max = 1.5;
    params.Nx = 400; params.Ny = 150;
    
    std::cout << "Params: ε=" << params.epsilon << ", α=" << params.alpha 
              << ", β₀=" << params.beta0 << std::endl;
    std::cout << "Domain: [" << params.Lx_min << "," << params.Lx_max 
              << "] × [" << params.Ly_min << "," << params.Ly_max << "]" << std::endl;
    std::cout << "λ₀=30.0, 100.0" << std::endl;
    std::cout << "Data saved in fig_lambda0_30.dat, fig_lambda0_100.dat" << std::endl;
    
    PhaseFieldModel model(params);
    model.run_parameter_effect_test_lambda0();
    std::cout << "λ₀ effect test completed, results saved successfully!" << std::endl;
}

void run_alpha_parameter_study() {
    std::cout << "\n=== α Effect ===" << std::endl;
    PhaseFieldModel::Parameters params;
    params.dx = 0.01; params.dy = 0.01;
    params.epsilon = 5 * params.dx / (4.0 * sqrt(2.0) * atanh(0.9));
    params.lambda0 = 30;
    params.beta0 = 0.5;
    params.total_time_steps = 40000;
    params.Lx_min = 0.0; params.Lx_max = 4.0;
    params.Ly_min = 0.0; params.Ly_max = 1.5;
    params.Nx = 400; params.Ny = 150;
    
    std::cout << "Params: ε=" << params.epsilon << ", λ₀=" << params.lambda0 
              << ", β₀=" << params.beta0 << std::endl;
     std::cout << "Domain: [" << params.Lx_min << "," << params.Lx_max 
              << "] × [" << params.Ly_min << "," << params.Ly_max << "]" << std::endl;
    std::cout << "α=0.05, 0.5, 1.0" << std::endl;

    
    PhaseFieldModel model(params);
    model.run_parameter_effect_test_alpha();
    std::cout << "α effect test completed, results saved successfully!" << std::endl;
}




int main() {
    try {
        auto total_start = std::chrono::high_resolution_clock::now();
        
        std::cout << "Program Started" << std::endl;
        
        PhaseFieldModel::Parameters temp_params;
        PhaseFieldModel temp_model(temp_params);


        std::cout << "========================================" << std::endl;
        std::cout << "      Phase-Field Calculating" << std::endl;
        std::cout << "========================================" << std::endl;

        run_convergence_analysis();

        run_lambda0_parameter_study();
    
        run_alpha_parameter_study();

        
        auto total_end = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
        
        std::cout << "\nStatistic:" << std::endl;
        std::cout << "Total Running Time: " << total_duration.count() << " sec" << std::endl;
        std::cout << "\n========================================" << std::endl;
        std::cout << "      All tasks completed!" << std::endl;
        std::cout << "========================================" << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "An error occurred: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cout << "An unexpected error occurred!" << std::endl;
        return 1;
    }

    
    
    return 0;
}