#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <sstream>
#include <time.h>
#include <omp.h>
#include "operations_omp.h"
#include "GradDescentOptim_omp.h"
#include "GCNLayer_omp.h"
#include "SoftmaxLayer_omp.h"
#include "GCNNetwork_omp.h"


using namespace std;


vector<double> cross_entropy(const vector<vector<double>>& pred, const vector<vector<double>>& labels) 
{
    vector<double> loss;
    
    int i,j;
    for (i = 0; i < pred.size(); i++) 
    {
        // Find the index of the maximum value in the ground truth label
        int maxIdx = 0;
        double maxVal = labels[i][0];
        for (j = 1; j < labels[i].size(); j++) 
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


void writeEmbedData(vector<vector<double>> data, string filepath)
{
    int rows = data.size();
    int cols = data[0].size();

    ofstream outFile(filepath);

    if (!outFile.is_open()) 
    {
        cerr << "Error opening file " << filepath << endl;
        return;
    }


    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j < cols; ++j) 
        {
            if(j!=cols-1)
            {
                outFile<<data[i][j]<<",";
            }
            else
            {
                outFile<<data[i][j];
            }
            
        }
        outFile<<endl;
        
    }
    cout<<"Successfully wrote the data as a binary file"<<endl;
    outFile.close();
}

void writeModelEval(vector<double> accs, vector<double> train_losses, vector<double> test_losses, string filepath)
{
    int n = accs.size();

    ofstream outFile(filepath);

    if (!outFile.is_open()) 
    {
        cerr << "Error opening file " << filepath << endl;
        return;
    }

    outFile<<"Iterations,Train loss,Test loss,Accuracy"<<endl;
    for(int i=0;i<n;i++)
    {
        outFile<<i+1<<","<<train_losses[i]<<","<<test_losses[i]<<","<<accs[i]<<endl;
    }
    cout<<"Successfully wrote the model evaluation parameters to a csv file !!"<<endl<<endl;
    outFile.close();
}


void writeTimeTaken(vector<int> num_threads, vector<double> time_taken, string filepath)
{
    int n = num_threads.size();

    ofstream outFile(filepath);

    if (!outFile.is_open()) 
    {
        cerr << "Error opening file " << filepath << endl;
        return;
    }

    outFile<<"Number of threads,Time taken"<<endl;

    for(int i=0;i<n;i++)
    {
        outFile<<num_threads[i]<<","<<time_taken[i]<<endl;
    }
    cout<<"Successfully wrote the time taken to a csv file !!"<<endl<<endl;
    outFile.close();
}




double GCN_train(int n,vector<vector<double>> A, int n_classes, vector<vector<double>> labels, int num_hidden_layers, vector<int> hidden_layer_sizes,int thread)
{
    omp_set_num_threads(thread);
    double timer = omp_get_wtime();
    cout<<endl<<"Training the model for "<<thread<<" number of threads and "<<num_hidden_layers<<" number of hidden layers"<<endl;
    int i,j,epoch;
    double sum;

    vector<vector<double>> A_mod = matrix_sum(A,identity_matrix(n));

    vector<vector<double>> D(n, vector<double>(n,0));
    

    for( i=0;i<n;i++)
    {
        sum = 0;
        for( j=0;j<A[0].size();j++)
        {
            sum+=A_mod[i][j];
        }
        D[i][i] = sum;
    }

    vector<vector<double>> D_mod = square_root_inverse_diag_matrix(D);
    
    vector<vector<double>> A_hat = matrix_product(matrix_product(D_mod,A_mod),D_mod);

    vector<vector<double>> X = identity_matrix(n);
    

    GCNNetwork gcn_model(n,n_classes,2,{16,2},activation_function,100);


    vector<int> train_nodes = {0,1,8};
    vector<int> test_nodes;

    for (int i = 0; i < labels.size(); ++i) 
    {
        if (find(train_nodes.begin(), train_nodes.end(), i) == train_nodes.end()) 
        {
            test_nodes.push_back(i);
        }
    }
    GradDescentOptim opt2(0.02,0.025);

    vector<vector<vector<double>>> embeds;
    vector<double> accs;
    vector<double> train_losses;
    vector<double> test_losses;

    double loss_min = pow(10,6);
    int es_iters = 0;
    int es_steps = 50;

    int max_epochs = 7000;

    double loss_train, loss_test;
    vector<double> loss;

    vector<bool> acc_test;


    for( epoch =0;epoch<max_epochs;epoch++)
    {
        loss_train=0;
        loss_test=0;
        vector<vector<double>> y_pred = gcn_model.forward(A_hat,X);

        opt2.operator()(y_pred,labels,train_nodes);

        gcn_model.sm_out->backward(opt2,true);

        for( i = gcn_model.gcn_layers.size()-1;i>=0;i--)
        {
            gcn_model.gcn_layers[i]->backward(opt2,true);
        }

        embeds.push_back(gcn_model.embedding(A_hat,X));

        vector<bool> acc_test;
        for (size_t i = 0; i < labels.size(); ++i) {
            if (find(train_nodes.begin(), train_nodes.end(), i) == train_nodes.end()) 
            {
                int argmax_y_pred = distance(y_pred[i].begin(), max_element(y_pred[i].begin(), y_pred[i].end()));
                int argmax_labels = distance(labels[i].begin(), max_element(labels[i].begin(), labels[i].end()));
                acc_test.push_back(argmax_y_pred == argmax_labels);
            }
        }
        double acc_mean = accumulate(acc_test.begin(), acc_test.end(), 0.0) / acc_test.size();
        accs.push_back(acc_mean);

        loss = cross_entropy(y_pred, labels);
        for( int i : train_nodes)
        {
            loss_train+=loss[i]/train_nodes.size();
        }
        for( int i : test_nodes)
        {
            loss_test+=loss[i]/test_nodes.size();
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

    timer = omp_get_wtime()-timer;
    
    cout<<endl<<"Finished training and classifying the data !!"<<endl<<endl;

    cout<<endl<<"Accuracy obtained : "<<accs.back()<<endl;

    cout<<"Total time taken for the program : "<<timer<<endl<<endl;

    string filepath1 = "SynData Solutions/omp/GCN_embed_omp_"+to_string(thread)+".csv";
    writeEmbedData(embeds.back(),filepath1);
    

    string filepath2 = "SynData Solutions/omp/GCN_omp_"+to_string(thread)+"_eval_params.csv";
    writeModelEval(accs,train_losses,test_losses,filepath2);


    
    return timer;
}




int main()
{

    
    
    string filename_1 = "Zacharys_Karate_club.csv";
    vector<vector<double>> A = readCSV(filename_1);
    int n = A.size();

    string filename_2 = "Communities.csv";
    vector<vector<double>> communities = readCSV(filename_2);
    int n_classes = communities.size();

    vector<int> num_threads = {1,2,4,6,8,10,12};

    vector<vector<double>> labels(n,vector<double>(communities.size(),0));

    for(int i=0;i<communities.size();i++)
    {
        for(int j=0;j<communities[i].size();j++)
        {
            labels[communities[i][j]][i] = 1;
        }
    }

    int num_hidden_layers = 2;
    vector<int> hidden_layer_sizes = {16,2};

    vector<double> time_taken;

    for(int i=0;i<num_threads.size();i++)
    {
        time_taken.push_back(GCN_train(n,A,n_classes,labels,num_hidden_layers,hidden_layer_sizes,num_threads[i]));
    }

    
    string filepath = "SynData Solutions/omp/GCN_omp_time_taken_"+to_string(num_hidden_layers)+".csv";
    writeTimeTaken(num_threads,time_taken,filepath);
    


    return 0;
    
}