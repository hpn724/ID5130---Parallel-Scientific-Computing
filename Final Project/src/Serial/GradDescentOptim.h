#ifndef GRADDESCENTOPTIM_H
#define GRADDESCENTOPTIM_H

#include <vector>

class GradDescentOptim {
private:
    double lr;
    double wd;
    std::vector<std::vector<double>> y_pred;
    std::vector<std::vector<double>> y_true;
    mutable std::vector<std::vector<double>> out;
    int bs;
    std::vector<int> train_nodes;

public:
    GradDescentOptim(double lr, double wd);

    void operator()(const std::vector<std::vector<double>>& y_pred, const std::vector<std::vector<double>>& y_true, const std::vector<int>& train_nodes = {});

    double getlr() const;
    double getwd() const;
    std::vector<std::vector<double>> getYPred() const;
    std::vector<std::vector<double>> getYTrue() const;
    std::vector<int> getTrainNodes() const;
    std::vector<std::vector<double>>& getOut() const;
    void setOut(const std::vector<std::vector<double>>& y) const;
    int getBatchSize() const;
};

#endif
