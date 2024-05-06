#ifndef SOFTMAXLAYER_H
#define SOFTMAXLAYER_H

#include <string>
#include <vector>
#include <functional>
#include "GradDescentOptim.h" // Include necessary header files

class SoftmaxLayer {
private:
    int n_inputs;
    int n_outputs;
    std::vector<std::vector<double>> W;
    std::vector<std::vector<double>> b;
    std::string name;
    std::vector<std::vector<double>> _X;

public:
    SoftmaxLayer(int n_inputs, int n_outputs, std::string name = "");

    std::string repr() const;

    std::vector<std::vector<double>> shift(const std::vector<std::vector<double>>& proj) const;

    std::vector<std::vector<double>> forward(const std::vector<std::vector<double>>& X, std::vector<std::vector<double>> W = {}, std::vector<std::vector<double>> b = {});

    std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> backward(const GradDescentOptim& optim, bool update = true);

    std::vector<std::vector<double>> getAttr(std::string argname);
};

#endif // SOFTMAXLAYER_H
