vector<vector<double>> dW_approx = get_grads({gcn1_out},sm1,"W",labels,0.0001,opt.getwd());
    vector<vector<double>> db_approx = get_grads({gcn1_out},sm1,"b",labels,0.0001,opt.getwd());
    
    pair<vector<vector<double>>, vector<vector<double>>> dWdb = sm1.backward(opt,false);

    vector<vector<double>> dW = dWdb.first;
    vector<vector<double>> db = dWdb.second;

    cout<<"dw : "<<endl;
    show_matrix(dW);
    cout<<"db : "<<endl;
    show_matrix(db);