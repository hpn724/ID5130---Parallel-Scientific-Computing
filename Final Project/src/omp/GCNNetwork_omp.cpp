#include "GCNNetwork_omp.h"
#include "GCNLayer_omp.h"
#include "SoftmaxLayer_omp.h"

GCNNetwork::GCNNetwork(int n_inputs, int n_outputs, int n_layers, std::vector<int> hidden_sizes, std::function<double(double)> activation, int seed)
    : n_inputs(n_inputs), n_outputs(n_outputs), n_layers(n_layers), hidden_sizes(hidden_sizes), activation(activation) 
{
    srand(seed);
    
    gcn_layers.reserve(n_layers);
    
    // Input layer
    std::unique_ptr<GCNLayer> gcn_in = std::make_unique<GCNLayer>(n_inputs, hidden_sizes[0], activation, "in");
    gcn_layers.push_back(std::move(gcn_in));
    
    // Hidden layers
    for (int layer = 0; layer < n_layers; ++layer) 
    {
        std::unique_ptr<GCNLayer> gcn = std::make_unique<GCNLayer>(gcn_layers.back()->getW().size(), hidden_sizes[layer], activation, "h" + std::to_string(layer));
        gcn_layers.push_back(std::move(gcn));
    }
    
    // Output layer
    sm_out = std::make_unique<SoftmaxLayer>(hidden_sizes.back(), n_outputs, "sm");
}

std::string GCNNetwork::repr() const {
    std::string repr_str;
    for (const auto& gcn_layer : gcn_layers) 
    {
        repr_str += gcn_layer->repr() + "\n";
    }
    repr_str += sm_out->repr() + "\n";
    return repr_str;
}

std::vector<std::vector<double>> GCNNetwork::embedding(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& X) {
    std::vector<std::vector<double>> H = X;
    for (const auto& gcn_layer : gcn_layers) 
    {
        H = gcn_layer->forward(A, H);
    }
    return H;
}

std::vector<std::vector<double>> GCNNetwork::forward(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& X) {
    std::vector<std::vector<double>> H = embedding(A, X);
    std::vector<std::vector<double>> p = sm_out->forward(H);
    return p;
}
