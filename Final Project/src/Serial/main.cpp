#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <sstream>
#include "operations.h"
#include "GradDescentOptim.h"
#include "GCNLayer.h"
#include "SoftmaxLayer.h"
#include "GCNNetwork.h"

using namespace std;


vector<double> cross_entropy(const vector<vector<double>>& pred, const vector<vector<double>>& labels) {
    vector<double> loss;
    
    for (size_t i = 0; i < pred.size(); i++) 
    {
        // Find the index of the maximum value in the ground truth label
        int maxIdx = 0;
        double maxVal = labels[i][0];
        for (size_t j = 1; j < labels[i].size(); j++) 
        {
            if (labels[i][j] > maxVal) 
            {
                maxVal = labels[i][j];
                maxIdx = j;
            }
        }
        
        // Calculate the cross-entropy loss
        double logPred = -log(pred[i][maxIdx]);
        
        loss.push_back(logPred);
    }
    
    return loss;
}


double activation_function(double x)
{
    return tanh(x);
}





vector<vector<double>> readCSV(const string& filename) 
{
    vector<vector<double>> matrix;
    ifstream file(filename);

    if (!file.is_open()) 
    {
        cerr << "Error opening file " << filename << endl;
        return matrix; // Return an empty matrix
    }

    string line;
    while (getline(file, line)) 
    {
        vector<double> row;
        stringstream ss(line);
        string cell;
        double value;

        while (getline(ss, cell, ',')) 
        {
            value = std::stod(cell); 
            row.push_back(value);
        }
        matrix.push_back(row);
    }

    file.close();
    return matrix;
}

int main()
{
    string filename_1 = "../../Zacharys_Karate_club.csv";
    vector<vector<double>> A = readCSV(filename_1);
    int n = A.size();

    string filename_2 = "../../Communities.csv";
    vector<vector<double>> communities = readCSV(filename_2);
    int n_classes = communities.size();

    vector<vector<double>> labels(n,vector<double>(communities.size(),0));

    for(int i=0;i<communities.size();i++)
    {
        for(int j=0;j<communities[i].size();j++)
        {
            labels[communities[i][j]][i] = 1;
        }
    }

    vector<vector<double>> A_mod = matrix_sum(A,identity_matrix(n));

    vector<vector<double>> D(n, vector<double>(n,0));
    
    for(int i=0;i<n;i++)
    {
        double sum = 0;
        for(int j=0;j<A[0].size();j++)
        {
            sum+=A_mod[i][j];
        }
        D[i][i] = sum;
    }

    vector<vector<double>> D_mod = square_root_inverse_diag_matrix(D);
    
    vector<vector<double>> A_hat = matrix_product(matrix_product(D_mod,A_mod),D_mod);

    vector<vector<double>> X = identity_matrix(n);

    GCNLayer gcn1(n,2,activation_function,"1");
    SoftmaxLayer sm1(2,n_classes,"SM");
    GradDescentOptim opt(0,1);
    
    GCNNetwork gcn_model(n,n_classes,2,{16,2},activation_function,100);

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /*                                      Training the model                                         */
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<int> train_nodes = {0,1,8};
    vector<int> test_nodes;
    for (int i = 0; i < labels.size(); ++i) 
    {
        if (find(train_nodes.begin(), train_nodes.end(), i) == train_nodes.end()) 
        {
            test_nodes.push_back(i);
        }
    }
    GradDescentOptim opt2(0.002,0.0025);

    vector<vector<vector<double>>> embeds;
    vector<double> accs;
    vector<double> train_losses;
    vector<double> test_losses;

    double loss_min = pow(10,6);
    int es_iters = 0;
    int es_steps = 50;

    int max_epochs = 20000;

    double loss_train, loss_test;
    vector<double> loss;

    vector<bool> acc_test;

    for(int epoch =0;epoch<max_epochs;epoch++)
    {
        loss_train=0;
        loss_test=0;
        vector<vector<double>> y_pred = gcn_model.forward(A_hat,X);

        opt2.operator()(y_pred,labels,train_nodes);

        gcn_model.sm_out->backward(opt2,true);

        for(int i = gcn_model.gcn_layers.size()-1;i>=0;i--)
        {
            gcn_model.gcn_layers[i]->backward(opt2,true);
        }

        vector<bool> acc_test;
        for (size_t i = 0; i < labels.size(); ++i) {
            if (find(train_nodes.begin(), train_nodes.end(), i) == train_nodes.end()) {
                int argmax_y_pred = distance(y_pred[i].begin(), max_element(y_pred[i].begin(), y_pred[i].end()));
                int argmax_labels = distance(labels[i].begin(), max_element(labels[i].begin(), labels[i].end()));
                acc_test.push_back(argmax_y_pred == argmax_labels);
            }
        }
        double acc_mean = accumulate(acc_test.begin(), acc_test.end(), 0.0) / acc_test.size();
        accs.push_back(acc_mean);

        loss = cross_entropy(y_pred, labels);
        for(int i=0;i<train_nodes.size();i++)
        {
            loss_train+=loss[train_nodes[i]]/train_nodes.size();
        }
        for(int i=0;i<test_nodes.size();i++)
        {
            loss_test+=loss[test_nodes[i]]/test_nodes.size();
        }

        train_losses.push_back(loss_train);
        test_losses.push_back(loss_test);

        if (loss_test<loss_min)
        {
            loss_min=loss_test;
            es_iters=0;
        }
        else
        {
            es_iters+=1;
        }

        
        if(epoch%100==0)
        {
            printf("Epoch : %d \t Train loss : %.3f \t Test loss : %.3f \n",epoch+1,loss_train,loss_test);
        }
        

    }


    return 0;
    
}