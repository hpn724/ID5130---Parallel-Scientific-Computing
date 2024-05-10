#ifndef GCNNETWORK_omp_H
#define GCNNETWORK_omp_H

#include <vector>
#include <memory>
#include <string>
#include <functional>
#include "GCNLayer_omp.h"
#include "SoftmaxLayer_omp.h"

class GCNNetwork {
private:
    int n_inputs;
    int n_outputs;
    int n_layers;
    std::vector<int> hidden_sizes;
    std::function<double(double)> activation;
    
public:
    std::vector<std::unique_ptr<GCNLayer>> gcn_layers;
    std::unique_ptr<SoftmaxLayer> sm_out;

    GCNNetwork(int n_inputs, int n_outputs, int n_layers, std::vector<int> hidden_sizes, std::function<double(double)> activation, int seed = 0);

    std::string repr() const;

    std::vector<std::vector<double>> embedding(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& X);

    std::vector<std::vector<double>> forward(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& X);
};

#endif // GCNNETWORK_H
