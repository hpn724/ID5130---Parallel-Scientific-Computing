#include "SoftmaxLayer_omp.h"
#include "operations_omp.h" 
#include <cmath>
#include <string>

using namespace std;

SoftmaxLayer::SoftmaxLayer(int n_inputs, int n_outputs, string name) 
    : n_inputs(n_inputs), n_outputs(n_outputs), name(name) 

{
    W = glorot_init(n_outputs, n_inputs);
    b = vector<vector<double>>(n_outputs, vector<double>(1, 0.0));
}

string SoftmaxLayer::repr() const 
{
    string temp = name.empty() ? "" : "_";
    return "Softmax: W" + temp + name + " (" + to_string(n_inputs) + ", " + to_string(n_outputs) + ")";
}

vector<vector<double>> SoftmaxLayer::shift(const vector<vector<double>>& proj) const 
{
    vector<vector<double>> shiftx = proj;
    vector<double> max_vals(proj.size(), 0.0);
    int i,j;
    #pragma omp parallel for collapse(2) shared(max_vals,proj) private(i,j)
    for (j = 0; j < proj[0].size(); ++j) 
    {
        for (i = 0; i < proj.size(); ++i) 
        {
            #pragma omp critical
            max_vals[i] = max(max_vals[i], proj[i][j]);
        }
    }
    
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < proj[0].size(); ++j) 
    {
        for (size_t i = 0; i < proj.size(); ++i) 
        {

            shiftx[i][j] =shiftx[i][j]-max_vals[i];
        }
    }

    
    vector<vector<double>> exps(proj.size(), vector<double>(proj[0].size(), 0.0));

    #pragma omp parallel for 
    for (size_t j = 0; j < proj[0].size(); ++j) 
    {
        double sum_exp = 0.0;
        for (size_t i = 0; i < proj.size(); ++i) 
        {
            exps[i][j] = exp(shiftx[i][j]);
            #pragma omp atomic
            sum_exp += exps[i][j];
        }
        for (size_t i = 0; i < proj.size(); ++i) 
        {
            #pragma omp atomic
            exps[i][j] /= sum_exp;
        }
    }
    return exps;
}

vector<vector<double>> SoftmaxLayer::forward(const vector<vector<double>>& X, vector<vector<double>> W, vector<vector<double>> b) 
{
    _X = transpose(X);
    vector<vector<double>> proj(n_outputs, vector<double>(_X[0].size(), 0.0));

    if(W.empty())
    {
        W = this->W;
    }
    if(b.empty())
    {
        b = this->b;
    }

    
    proj = matrix_product(W,_X);
    
    proj = matrix_sum(proj,b);

    return transpose(shift(proj));
}

pair<vector<vector<double>>, vector<vector<double>>> SoftmaxLayer::backward(const GradDescentOptim& optim, bool update) 
{
    // Build mask on loss
    vector<double> train_mask_1(optim.getYPred().size(), 0.0);
    #pragma omp parallel for
    for (int node : optim.getTrainNodes()) 
    {
        train_mask_1[node] = 1.0;
    }

    vector<vector<double>> train_mask = reshape(train_mask_1,train_mask_1.size(),1);
    

    // Derivative of loss w.r.t. activation (pre-softmax)
    vector<vector<double>> d1(optim.getYPred().size(), vector<double>(n_outputs, 0.0));
    
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < optim.getYPred().size(); ++i) 
    {
        for (int j = 0; j < n_outputs; ++j) 
        {
            d1[i][j] = optim.getYPred()[i][j] - optim.getYTrue()[i][j];
            d1[i][j] *= train_mask[i][0];
        }

    }
    
    vector<vector<double>> grad = matrix_product(d1,W);
    optim.setOut(grad);

    vector<vector<double>> dW = matrix_product(transpose(d1),transpose(_X));

    vector<vector<double>> db(n_outputs, vector<double>(1, 0.0));

    #pragma omp parallel for collapse(2)
    for(int i=0;i<dW.size();i++)
    {
        for(int j=0;j<dW[0].size();j++)
        {
            dW[i][j] = dW[i][j]/optim.getBatchSize();
        }
    }
    vector<vector<double>> d1t = transpose(d1);
    
    int i,j;
    double sum;
    for(i=0; i<db.size();i++)
    {
        sum=0;
        for(j=0;j<d1t[i].size();j++)
        {
            #pragma omp atomic
            sum+=d1t[i][j];
        }
        
        db[i][0] = sum/optim.getBatchSize();
    }
    

    if (update) 
    {
        #pragma omp parallel for 
        for (size_t i = 0; i < n_outputs; ++i) 
        {
            for (int j = 0; j < n_inputs; ++j) 
            {
                #pragma omp atomic
                W[i][j] -= (dW[i][j] + W[i][j] * optim.getwd() / optim.getBatchSize()) * optim.getlr();
            }
            #pragma omp atomic
            b[i][0] -= (db[i][0] + b[i][0] * optim.getwd() / optim.getBatchSize()) * optim.getlr();
        }
    }

    return {dW, db};
}

vector<vector<double>> SoftmaxLayer::getAttr(string argname) 
{
    if(argname == "W")
    {
        return W;
    }
    else
    {
        
        return b;
    }
}