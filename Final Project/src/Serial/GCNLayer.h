#ifndef GCNLAYER_H
#define GCNLAYER_H

#include <vector>
#include <string>
#include <functional>
#include "GradDescentOptim.h" // Include necessary headers

class GCNLayer {
private:
    int n_inputs;
    int n_outputs;
    std::vector<std::vector<double>> W;
    std::function<double(double)> activation;
    std::string name;

    std::vector<std::vector<double>> _A;
    std::vector<std::vector<double>> _X;
    std::vector<std::vector<double>> _H;
    std::vector<std::vector<double>> grad;

public:
    GCNLayer(int n_inputs, int n_outputs, std::function<double(double)> activation, std::string name = "");

    std::string repr() const;

    std::vector<std::vector<double>> forward(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& X, std::vector<std::vector<double>> W = {});

    std::vector<std::vector<double>> backward(const GradDescentOptim& optim, bool update = true);

    std::vector<std::vector<double>> getW();
};

#endif
