#include "GradDescentOptim.h"

GradDescentOptim::GradDescentOptim(double lr, double wd) : lr(lr), wd(wd), bs(0) {}

void GradDescentOptim::operator()(const std::vector<std::vector<double>>& y_pred, const std::vector<std::vector<double>>& y_true, const std::vector<int>& train_nodes) {
    this->y_pred = y_pred;
    this->y_true = y_true;

    if (train_nodes.empty()) 
    {
        this->train_nodes.resize(y_pred.size());
        for (int i = 0; i < y_pred.size(); ++i) 
        {
            this->train_nodes[i] = i;
        }
    } else {
        this->train_nodes = train_nodes;
    }

    this->bs = this->train_nodes.size();
}

double GradDescentOptim::getlr() const 
{
    return lr;
}

double GradDescentOptim::getwd() const 
{
    return wd;
}

std::vector<std::vector<double>> GradDescentOptim::getYPred() const 
{
    return y_pred;
}

std::vector<std::vector<double>> GradDescentOptim::getYTrue() const 
{
    return y_true;
}

std::vector<int> GradDescentOptim::getTrainNodes() const 
{
    return train_nodes;
}

std::vector<std::vector<double>>& GradDescentOptim::getOut() const 
{
    return out;
}

void GradDescentOptim::setOut(const std::vector<std::vector<double>>& y) const 
{
    this->out = y;
}

int GradDescentOptim::getBatchSize() const 
{
    return this->bs;
}
