#include "GCNLayer_omp.h"
#include "GradDescentOptim_omp.h" 
#include "operations_omp.h"
#include<random>
#include<cmath>
#include<omp.h>

using namespace std;


GCNLayer::GCNLayer(int n_inputs, int n_outputs, std::function<double(double)> activation, std::string name)
    : n_inputs(n_inputs), n_outputs(n_outputs), activation(activation), name(name) 
    {
    W = glorot_init(n_outputs, n_inputs); // Assuming glorot_init is defined elsewhere
    }

std::string GCNLayer::repr() const 
{
    std::string temp;
    if (name.empty()) 
    {
        temp = "";
    } else 
    {
        temp = "_";
    }
    return "GCN: W" + temp + name + " (" + std::to_string(n_inputs) + ", " + std::to_string(n_outputs) + ")";
}

std::vector<std::vector<double>> GCNLayer::forward(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& X, std::vector<std::vector<double>> W) 
{
    _A = A;
    _X = transpose(matrix_product(A,X));


    if (W.empty()) 
    {
        W = this->W;
    }
    
    vector<vector<double>> H = matrix_product(W,_X);
    
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < H.size(); ++i) 
    {
        for (size_t j = 0; j < H[0].size(); ++j) 
        {
            if (activation != nullptr) 
            {
                H[i][j] = activation(H[i][j]);
            }
        }
    }

    this->_H = H;

    return transpose(this->_H);
}

std::vector<std::vector<double>> GCNLayer::backward(const GradDescentOptim& optim, bool update) 
{
    vector<vector<double>> dtanh(_H.size(), vector<double>(_H[0].size()));

    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < _H.size(); ++i) 
    {
        for (size_t j = 0; j < _H[0].size(); ++j) 
        {
            dtanh[i][j] = 1 - pow(_H[i][j], 2);
        }
    }


    vector<vector<double>> d2(optim.getOut().size(), vector<double>(optim.getOut()[0].size()));

    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < optim.getOut().size(); ++i) 
    {
        for (size_t j = 0; j < optim.getOut()[0].size(); ++j) 
        {
            d2[i][j] = optim.getOut()[i][j] * dtanh[j][i];
        }
    }
    

    vector<vector<double>> grad = matrix_product(matrix_product(_A,d2),W);
    
    optim.setOut(grad);
    

    vector<vector<double>> dW = matrix_product(transpose(d2),transpose(_X));
    #pragma omp parallel for collapse(2)
    for(int i=0;i<dW.size();i++)
    {
        for(int j=0;j<dW[0].size();j++)
        {
            dW[i][j] = dW[i][j]/optim.getBatchSize();
        }
    }
    
    vector<vector<double>> dW_wd(W.size(), vector<double>(W[0].size()));

    #pragma omp parallel for collapse(2)
    for(int i=0;i<W.size();i++)
    {
        for(int j=0;j<W[0].size();j++)
        {
            dW_wd[i][j] = W[i][j]*optim.getwd()/optim.getBatchSize();
        }
    }

    if (update) 
    {
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < W.size(); ++i) 
        {
            for (size_t j = 0; j < W[0].size(); ++j) 
            {
                #pragma omp atomic
                W[i][j] -= (dW[i][j] + dW_wd[i][j]) * optim.getlr();
            }
        }
    }

    return dW;
}

std::vector<std::vector<double>> GCNLayer::getW() 
{
    return W;
}
