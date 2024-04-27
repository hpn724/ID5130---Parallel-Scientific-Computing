import matplotlib.pyplot as plt
import pandas as pd



###########################################################################
###                           Question 1                                ###
###########################################################################

############################    Part - a    ###############################


wave_1D_serial_upwind = pd.read_csv('Solution/Q1a/1D_Wave_serial_upwind.csv',index_col=0,header=0)
wave_1D_serial_QUICK = pd.read_csv('Solution/Q1a/1D_Wave_serial_QUICK.csv',index_col=0,header=0)
wave_1D_serial_anal = pd.read_csv('Solution/Q1a/1D_Wave_serial_analytical.csv',index_col=0,header=0)

x1 = list(wave_1D_serial_upwind.columns.astype(float))

wave_1D_serial_upwind_t0 = wave_1D_serial_upwind.loc[0.0][:].tolist()
wave_1D_serial_upwind_t1 = wave_1D_serial_upwind.loc[0.5][:].tolist()
wave_1D_serial_upwind_t2 = wave_1D_serial_upwind.loc[1.0][:].tolist()

wave_1D_serial_QUICK_t0 = wave_1D_serial_QUICK.loc[0.0][:].tolist()
wave_1D_serial_QUICK_t1 = wave_1D_serial_QUICK.loc[0.5][:].tolist()
wave_1D_serial_QUICK_t2 = wave_1D_serial_QUICK.loc[1.0][:].tolist()

wave_1D_serial_anal_t0 = wave_1D_serial_anal.loc[0.0][:].tolist()
wave_1D_serial_anal_t1 = wave_1D_serial_anal.loc[0.5][:].tolist()
wave_1D_serial_anal_t2 = wave_1D_serial_anal.loc[1.0][:].tolist()



oneD_serial_t0_fig = plt.figure()
plt.plot(x1,wave_1D_serial_anal_t0,'r',label='analytical',linestyle='dotted')
plt.plot(x1,wave_1D_serial_upwind_t0,'b-',label='upwind')
plt.plot(x1,wave_1D_serial_QUICK_t0,'g-',label='QUICK')
plt.grid()
plt.legend()
plt.title('One-Dimensional Wave solution for t=0')
plt.savefig("Solution/Q1a/1D_wave_serial_t0.png")

oneD_serial_t1_fig = plt.figure()
plt.plot(x1,wave_1D_serial_anal_t1,'r',label='analytical',linestyle='dotted')
plt.plot(x1,wave_1D_serial_upwind_t1,'b-',label='upwind')
plt.plot(x1,wave_1D_serial_QUICK_t1,'g-',label='QUICK')
plt.grid()
plt.legend()
plt.title('One-Dimensional Wave solution for t=0.5')
plt.savefig("Solution/Q1a/1D_wave_serial_t1.png")

oneD_serial_t2_fig = plt.figure()
plt.plot(x1,wave_1D_serial_anal_t2,'r',label='analytical',linestyle='dotted')
plt.plot(x1,wave_1D_serial_upwind_t2,'b-',label='upwind')
plt.plot(x1,wave_1D_serial_QUICK_t2,'g-',label='QUICK')
plt.grid()
plt.legend()
plt.title('One-Dimensional Wave solution for t=1')
plt.savefig("Solution/Q1a/1D_wave_serial_t2.png")



############################    Part - b    ###############################



num_procs_Q1 = [2,4]
for procs in num_procs_Q1:
    
    files = ['upwind','QUICK','analytical']
    filename_upwind = '1D_Wave_MPI_'+str(procs)+'_upwind.csv'
    filename_QUICK = '1D_Wave_MPI_'+str(procs)+'_QUICK.csv'
    filename_analytical = '1D_Wave_MPI_'+str(procs)+'_analytical.csv'
    wave_1D_MPI_upwind = pd.read_csv('Solution/Q1b/'+filename_upwind,index_col=0,header=0)
    wave_1D_MPI_QUICK = pd.read_csv('Solution/Q1b/'+filename_QUICK,index_col=0,header=0)
    wave_1D_MPI_anal = pd.read_csv('Solution/Q1b/'+filename_analytical,index_col=0,header=0)

    x1 = list(wave_1D_MPI_upwind.columns.astype(float))

    wave_1D_MPI_upwind_t0 = wave_1D_MPI_upwind.loc[0.0][:].tolist()
    wave_1D_MPI_upwind_t1 = wave_1D_MPI_upwind.loc[0.5][:].tolist()
    wave_1D_MPI_upwind_t2 = wave_1D_MPI_upwind.loc[1.0][:].tolist()

    wave_1D_MPI_QUICK_t0 = wave_1D_MPI_QUICK.loc[0.0][:].tolist()
    wave_1D_MPI_QUICK_t1 = wave_1D_MPI_QUICK.loc[0.5][:].tolist()
    wave_1D_MPI_QUICK_t2 = wave_1D_MPI_QUICK.loc[1.0][:].tolist()

    wave_1D_MPI_anal_t0 = wave_1D_MPI_anal.loc[0.0][:].tolist()
    wave_1D_MPI_anal_t1 = wave_1D_MPI_anal.loc[0.5][:].tolist()
    wave_1D_MPI_anal_t2 = wave_1D_MPI_anal.loc[1.0][:].tolist()



    oneD_MPI_t0_fig = plt.figure()
    plt.plot(x1,wave_1D_MPI_anal_t0,'r',label='analytical',linestyle='dotted')
    plt.plot(x1,wave_1D_MPI_upwind_t0,'b-',label='upwind')
    plt.plot(x1,wave_1D_MPI_QUICK_t0,'g-',label='QUICK')
    plt.grid()
    plt.legend()
    plt.title('One-Dimensional Wave solution for t=0, procs = '+str(procs))
    plt.savefig("Solution/Q1b/1D_wave_MPI_"+str(procs)+"_t0.png")

    oneD_MPI_t1_fig = plt.figure()
    plt.plot(x1,wave_1D_MPI_anal_t1,'r',label='analytical',linestyle='dotted')
    plt.plot(x1,wave_1D_MPI_upwind_t1,'b-',label='upwind')
    plt.plot(x1,wave_1D_MPI_QUICK_t1,'g-',label='QUICK')
    plt.grid()
    plt.legend()
    plt.title('One-Dimensional Wave solution for t=0.5, procs = '+str(procs))
    plt.savefig("Solution/Q1b/1D_wave_MPI_"+str(procs)+"_t1.png")

    oneD_MPI_t2_fig = plt.figure()
    plt.plot(x1,wave_1D_MPI_anal_t2,'r',label='analytical',linestyle='dotted')
    plt.plot(x1,wave_1D_MPI_upwind_t2,'b-',label='upwind')
    plt.plot(x1,wave_1D_MPI_QUICK_t2,'g-',label='QUICK')
    plt.grid()
    plt.legend()
    plt.title('One-Dimensional Wave solution for t=1, procs = '+str(procs))
    plt.savefig("Solution/Q1b/1D_wave_MPI_"+str(procs)+"_t2.png")
    

###########################################################################
###                             Question 2                              ###
###########################################################################


############################    Part - a    ###############################

poisson_jacobi_serial = pd.read_csv('Solution/Q2a/Poisson_Jacobi_serial.csv',index_col=0,header=0)

x_poisson_serial = poisson_jacobi_serial.index.astype(float).tolist()
y_poisson_serial = poisson_jacobi_serial.columns.astype(float).tolist()

phi_y_poisson_jacobi_serial = poisson_jacobi_serial.loc[0.0][:].tolist()
phi_x_poisson_jacobi_serial = poisson_jacobi_serial['0'].tolist()

poisson_jacobi_serial_x_fig = plt.figure()
plt.plot(x_poisson_serial,phi_x_poisson_jacobi_serial,label='phi vs x, y=0')
plt.xlabel('x coordinates')
plt.ylabel('phi value')
plt.legend()
plt.title('phi vs x  , Poisson Jacobi iteration serial program')
plt.grid()
plt.savefig("Solution/Q2a/phi_vs_x_Poisson_Jacobi_serial.png")

poisson_jacobi_serial_y_fig = plt.figure()
plt.plot(y_poisson_serial,phi_y_poisson_jacobi_serial,label='phi vs y, x=0')
plt.xlabel('y coordinates')
plt.ylabel('phi value')
plt.legend()
plt.title('phi vs y  , Poisson Jacobi iteration serial program')
plt.grid()
plt.savefig("Solution/Q2a/phi_vs_y_Poisson_Jacobi_serial.png")



############################    Part - b    ###############################

poisson_jacobi_MPI_2_01 = pd.read_csv('Solution/Q2b/Poisson_Jacobi_MPI_2_0.010.csv',index_col=0,header=0)
poisson_jacobi_MPI_4_01 = pd.read_csv('Solution/Q2b/Poisson_Jacobi_MPI_4_0.010.csv',index_col=0,header=0)
poisson_jacobi_MPI_6_01 = pd.read_csv('Solution/Q2b/Poisson_Jacobi_MPI_6_0.010.csv',index_col=0,header=0)
poisson_jacobi_serial_01 = pd.read_csv('Solution/Q2b/Poisson_Jacobi_serial0.010.csv',index_col=0,header=0)

x_poisson_MPI_2_01 = poisson_jacobi_MPI_2_01.index.astype(float).tolist()
y_poisson_MPI_2_01 = poisson_jacobi_MPI_2_01.columns.astype(float).tolist()

x_poisson_MPI_4_01 = poisson_jacobi_MPI_4_01.index.astype(float).tolist()
y_poisson_MPI_4_01 = poisson_jacobi_MPI_4_01.columns.astype(float).tolist()

x_poisson_MPI_6_01 = poisson_jacobi_MPI_6_01.index.astype(float).tolist()
y_poisson_MPI_6_01 = poisson_jacobi_MPI_6_01.columns.astype(float).tolist()

x_poisson_serial_01 = poisson_jacobi_serial_01.index.astype(float).tolist()
y_poisson_serial_01 = poisson_jacobi_serial_01.columns.astype(float).tolist()

phi_y_poisson_jacobi_MPI_2_01 = poisson_jacobi_MPI_2_01.loc[0.0][:].tolist()
phi_x_poisson_jacobi_MPI_2_01 = poisson_jacobi_MPI_2_01['0'].tolist()

phi_y_poisson_jacobi_MPI_4_01 = poisson_jacobi_MPI_4_01.loc[0.0][:].tolist()
phi_x_poisson_jacobi_MPI_4_01 = poisson_jacobi_MPI_4_01['0'].tolist()

phi_y_poisson_jacobi_MPI_6_01 = poisson_jacobi_MPI_6_01.loc[0.0][:].tolist()
phi_x_poisson_jacobi_MPI_6_01 = poisson_jacobi_MPI_6_01['0'].tolist()

phi_y_poisson_jacobi_serial_01 = poisson_jacobi_serial_01.loc[0.0][:].tolist()
phi_x_poisson_jacobi_serial_01 = poisson_jacobi_serial_01['0'].tolist()

poisson_jacobi_MPI_01_x_fig = plt.figure()
plt.plot(x_poisson_MPI_2_01,phi_x_poisson_jacobi_MPI_2_01,label='processors = 2')
plt.plot(x_poisson_MPI_4_01,phi_x_poisson_jacobi_MPI_4_01,label='processors = 4')
plt.plot(x_poisson_MPI_6_01,phi_x_poisson_jacobi_MPI_6_01,label='processors = 6')
plt.plot(x_poisson_serial_01,phi_x_poisson_jacobi_serial_01,label='serial',linestyle='dotted')
plt.xlabel('x coordinates')
plt.ylabel('phi value')
plt.legend()
plt.title('phi vs x, y = 0  , Poisson Jacobi iteration MPI program grid size = 0.01')
plt.grid()
plt.savefig("Solution/Q2b/phi_vs_x_Poisson_Jacobi_MPI_01.png")

poisson_jacobi_MPI_01_y_fig = plt.figure()
plt.plot(y_poisson_MPI_2_01,phi_y_poisson_jacobi_MPI_2_01,label='processors = 2')
plt.plot(y_poisson_MPI_4_01,phi_y_poisson_jacobi_MPI_4_01,label='processors = 4')
plt.plot(y_poisson_MPI_6_01,phi_y_poisson_jacobi_MPI_6_01,label='processors = 6')
plt.plot(y_poisson_serial_01,phi_y_poisson_jacobi_serial_01,label='serial',linestyle='dotted')
plt.xlabel('y coordinates')
plt.ylabel('phi value')
plt.legend()
plt.title('phi vs y, x = 0  , Poisson Jacobi iteration MPI program grid size = 0.01')
plt.grid()
plt.savefig("Solution/Q2b/phi_vs_y_Poisson_Jacobi_MPI_01.png")


poisson_jacobi_MPI_2_iter_01 = pd.read_csv('Solution/Q2b/Poisson_Jacobi_error_MPI_2_0.010.csv',header = 0)
poisson_jacobi_MPI_4_iter_01 = pd.read_csv('Solution/Q2b/Poisson_Jacobi_error_MPI_4_0.010.csv',header = 0)
poisson_jacobi_MPI_6_iter_01 = pd.read_csv('Solution/Q2b/Poisson_Jacobi_error_MPI_6_0.010.csv',header = 0)

poisson_jacobi_MPI_01_error = plt.figure()
plt.plot(poisson_jacobi_MPI_2_iter_01['Iteration'],poisson_jacobi_MPI_2_iter_01['Error'], label = 'processor = 2')
plt.plot(poisson_jacobi_MPI_4_iter_01['Iteration'],poisson_jacobi_MPI_4_iter_01['Error'], label = 'processor = 4')
plt.plot(poisson_jacobi_MPI_6_iter_01['Iteration'],poisson_jacobi_MPI_6_iter_01['Error'], label = 'processor = 6')
plt.xlabel('iterations')
plt.ylabel('error')
plt.legend()
plt.title('error vs iteration plot  , Poisson Jacobi iteration MPI program grid size = 0.01')
plt.grid()
plt.savefig("Solution/Q2b/error_Poisson_Jacobi_MPI_01.png")






############################    Part - c    ###############################

Poisson_GS_red_black_serial = pd.read_csv('Solution/Q2c/Poisson_GS_red_black_serial0.100.csv',index_col=0,header=0)

x_Poisson_GS_red_black_serial = Poisson_GS_red_black_serial.index.astype(float).tolist()
y_Poisson_GS_red_black_serial = Poisson_GS_red_black_serial.columns.astype(float).tolist()

phi_y_Poisson_GS_red_black_serial = Poisson_GS_red_black_serial.loc[0.0][:].tolist()
phi_x_Poisson_GS_red_black_serial = Poisson_GS_red_black_serial['0'].tolist()

Poisson_GS_red_black_serial_x_fig = plt.figure()
plt.plot(x_Poisson_GS_red_black_serial,phi_x_Poisson_GS_red_black_serial,label='phi vs x, y=0')
plt.xlabel('x coordinates')
plt.ylabel('phi value')
plt.legend()
plt.title('phi vs x  , Poisson Gauss Seidel Red Black iteration serial program')
plt.grid()
plt.savefig("Solution/Q2c/phi_vs_x_Poisson_GS_red_black_serial_1.png")


Poisson_GS_red_black_serial_y_fig = plt.figure()
plt.plot(y_Poisson_GS_red_black_serial,phi_y_Poisson_GS_red_black_serial,label='phi vs y, x=0')
plt.xlabel('y coordinates')
plt.ylabel('phi value')
plt.legend()
plt.title('phi vs y  , Poisson Gauss Seidel Red Black iteration serial program')
plt.grid()
plt.savefig("Solution/Q2c/phi_vs_y_Poisson_GS_red_black_serial_1.png")


poisson_GS_red_black_MPI_2_01 = pd.read_csv('Solution/Q2c/Poisson_GS_red_black_MPI_2_0.010.csv',index_col=0,header=0)
poisson_GS_red_black_MPI_4_01 = pd.read_csv('Solution/Q2c/Poisson_GS_red_black_MPI_4_0.010.csv',index_col=0,header=0)
poisson_GS_red_black_MPI_6_01 = pd.read_csv('Solution/Q2c/Poisson_GS_red_black_MPI_6_0.010.csv',index_col=0,header=0)
poisson_GS_red_black_serial_01 = pd.read_csv('Solution/Q2c/Poisson_GS_red_black_serial0.010.csv',index_col=0,header=0)

x_poisson_MPI_2_01 = poisson_GS_red_black_MPI_2_01.index.astype(float).tolist()
y_poisson_MPI_2_01 = poisson_GS_red_black_MPI_2_01.columns.astype(float).tolist()

x_poisson_MPI_4_01 = poisson_GS_red_black_MPI_4_01.index.astype(float).tolist()
y_poisson_MPI_4_01 = poisson_GS_red_black_MPI_4_01.columns.astype(float).tolist()

x_poisson_MPI_6_01 = poisson_GS_red_black_MPI_6_01.index.astype(float).tolist()
y_poisson_MPI_6_01 = poisson_GS_red_black_MPI_6_01.columns.astype(float).tolist()

x_poisson_serial_01 = poisson_GS_red_black_serial_01.index.astype(float).tolist()
y_poisson_serial_01 = poisson_GS_red_black_serial_01.columns.astype(float).tolist()

phi_y_poisson_GS_red_black_MPI_2_01 = poisson_GS_red_black_MPI_2_01.loc[0.0][:].tolist()
phi_x_poisson_GS_red_black_MPI_2_01 = poisson_GS_red_black_MPI_2_01['0'].tolist()

phi_y_poisson_GS_red_black_MPI_4_01 = poisson_GS_red_black_MPI_4_01.loc[0.0][:].tolist()
phi_x_poisson_GS_red_black_MPI_4_01 = poisson_GS_red_black_MPI_4_01['0'].tolist()

phi_y_poisson_GS_red_black_MPI_6_01 = poisson_GS_red_black_MPI_6_01.loc[0.0][:].tolist()
phi_x_poisson_GS_red_black_MPI_6_01 = poisson_GS_red_black_MPI_6_01['0'].tolist()

phi_y_poisson_GS_red_black_serial_01 = poisson_GS_red_black_serial_01.loc[0.0][:].tolist()
phi_x_poisson_GS_red_black_serial_01 = poisson_GS_red_black_serial_01['0'].tolist()

poisson_GS_red_black_MPI_01_x_fig = plt.figure()
plt.plot(x_poisson_MPI_2_01,phi_x_poisson_GS_red_black_MPI_2_01,label='processors = 2')
plt.plot(x_poisson_MPI_4_01,phi_x_poisson_GS_red_black_MPI_4_01,label='processors = 4')
plt.plot(x_poisson_MPI_6_01,phi_x_poisson_GS_red_black_MPI_6_01,label='processors = 6')
plt.plot(x_poisson_serial_01,phi_x_poisson_GS_red_black_serial_01,label='serial',linestyle='dotted')
plt.xlabel('x coordinates')
plt.ylabel('phi value')
plt.legend()
plt.title('phi vs x, y = 0  , Poisson GS_red_black iteration MPI program grid size = 0.01')
plt.grid()
plt.savefig("Solution/Q2c/phi_vs_x_Poisson_GS_red_black_MPI_01.png")

poisson_GS_red_black_MPI_01_y_fig = plt.figure()
plt.plot(y_poisson_MPI_2_01,phi_y_poisson_GS_red_black_MPI_2_01,label='processors = 2')
plt.plot(y_poisson_MPI_4_01,phi_y_poisson_GS_red_black_MPI_4_01,label='processors = 4')
plt.plot(y_poisson_MPI_6_01,phi_y_poisson_GS_red_black_MPI_6_01,label='processors = 6')
plt.plot(y_poisson_serial_01,phi_y_poisson_GS_red_black_serial_01,label='serial',linestyle='dotted')
plt.xlabel('y coordinates')
plt.ylabel('phi value')
plt.legend()
plt.title('phi vs y, x = 0  , Poisson GS_red_black iteration MPI program grid size = 0.01')
plt.grid()
plt.savefig("Solution/Q2c/phi_vs_y_Poisson_GS_red_black_MPI_01.png")


poisson_GS_red_black_MPI_2_iter_01 = pd.read_csv('Solution/Q2c/Poisson_GS_red_black_error_MPI_2_0.010.csv',header = 0)
poisson_GS_red_black_MPI_4_iter_01 = pd.read_csv('Solution/Q2c/Poisson_GS_red_black_error_MPI_4_0.010.csv',header = 0)
poisson_GS_red_black_MPI_6_iter_01 = pd.read_csv('Solution/Q2c/Poisson_GS_red_black_error_MPI_6_0.010.csv',header = 0)

poisson_GS_red_black_MPI_01_error = plt.figure()
plt.plot(poisson_GS_red_black_MPI_2_iter_01['Iteration'],poisson_GS_red_black_MPI_2_iter_01['Error'], label = 'processor = 2')
plt.plot(poisson_GS_red_black_MPI_4_iter_01['Iteration'],poisson_GS_red_black_MPI_4_iter_01['Error'], label = 'processor = 4')
plt.plot(poisson_GS_red_black_MPI_6_iter_01['Iteration'],poisson_GS_red_black_MPI_6_iter_01['Error'], label = 'processor = 6')
plt.xlabel('iterations')
plt.ylabel('error')
plt.legend()
plt.title('error vs iteration plot  , Poisson GS_red_black iteration MPI program grid size = 0.01')
plt.grid()
plt.savefig("Solution/Q2c/error_Poisson_GS_red_black_MPI_01.png")



############################    Part - d    ###############################

poisson_jacobi_MPI_005 = pd.read_csv('Solution/Q2d/Poisson_Jacobi_speedup.csv',header=None)
poisson_GS_red_black_MPI_005 = pd.read_csv('Solution/Q2d/Poisson_GS_red_black_speedup.csv',header=None)

poisson_jacobi_MPI_005.rename(columns={0:'Processor',1:'Speedup'},inplace=True)
poisson_GS_red_black_MPI_005.rename(columns={0:'Processor',1:'Speedup'},inplace=True)

speedup_fig = plt.figure()
plt.plot(poisson_jacobi_MPI_005['Processor'],poisson_jacobi_MPI_005['Speedup'],label='Jacobi')
plt.plot(poisson_GS_red_black_MPI_005['Processor'],poisson_GS_red_black_MPI_005['Speedup'],label='Gauss-Seidel Red black')
plt.xlabel('Number of processors')
plt.ylabel('Speedup')
plt.title('Speedup plot')
plt.grid()
plt.legend()
plt.savefig("Solution/Q2d/speedup_plot.png")