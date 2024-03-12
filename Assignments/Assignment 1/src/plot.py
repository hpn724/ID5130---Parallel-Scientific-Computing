import matplotlib.pyplot as plt
import pandas as pd

LU_data = pd.read_csv('src/Pade_scheme_LU_Decomp_plot_data.csv')

RD_data = pd.read_csv('src/Pade_scheme_RecursiveDoubling_data.csv')
RD_time_data=pd.read_csv('src/PadeScheme_RecursiveDoubling_Thread_Time.csv')

GS_serial_data=pd.read_csv('src/Poisson_GaussSeidel_serial_data.csv')
GS_Diag_data_g01_thr2=pd.read_csv('src/Poisson_GaussSeidel_DiagonalWave_data_grid_0.100000_threads_2.csv')
GS_BR_data_g01_thr2=pd.read_csv('src/Poisson_GaussSeidel_RedBlack_data_grid_0.100000_threads_2.csv')

GS_DiagWave_grid_time=pd.read_csv('src/Poisson_GaussSeidel_DiagonalWave_data_grid_time_data.csv')
GS_RedBlack_grid_time=pd.read_csv('src/Poisson_GaussSeidel_RedBlack_data_grid_time_data.csv')

GS_DiagWave_thread_time=pd.read_csv('src/Poisson_GaussSeidel_DiagonalWave_data_thread_time_data.csv')
GS_RedBlack_thread_time=pd.read_csv('src/Poisson_GaussSeidel_RedBlack_data_thread_time_data.csv')


plot1 = plt.figure()
plt.plot(LU_data['x'],LU_data['Actual value'],'o',label='5cos(x)')
plt.plot(LU_data['x'],LU_data['LU Decomposition values'],label='LU Decomp')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('LU Decomposition')
plt.grid()
fig1=plot1.get_figure()
fig1.savefig('Result_Images_and_Docs/Pade_LU_Decomp.png')

plot2 = plt.figure()
plt.plot(RD_data['x'],RD_data['Actual value'],'o',label='5cos(x)')
plt.plot(RD_data['x'],RD_data['Recursive Doubling'],label='Recursive Doubling')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.title('Pade scheme Recursive Doubling')
fig2=plot2.get_figure()
fig2.savefig('Result_Images_and_Docs/Pade_scheme_RecursiveDoubling.png')

plot3 = plt.figure()
plt.plot(RD_time_data['Thread Count'],RD_time_data['Time'])
plt.xlabel('No of threads')
plt.ylabel('Time taken')
plt.title('Pade scheme Recursive Doubling time taken per thread')
plt.grid()
fig3=plot3.get_figure()
fig3.savefig('Result_Images_and_Docs/Pade_scheme_RecursiveDoubling_time_taken.png')


plot4 = plt.figure()
plt.plot(GS_serial_data['x'],GS_serial_data['y values'],'o',label='2(2-x^2-y^2)')
plt.plot(GS_serial_data['x'],GS_serial_data['y solved'],label='Gauss Seidel Serial')
plt.xlabel('x')
plt.ylabel('f(x,y)')
plt.legend()
plt.grid()
plt.title('Gauss Seidel Serial data for y=0.5')
fig4=plot4.get_figure()
fig4.savefig('Result_Images_and_Docs/Gauss_Seidel_Serial.png')


plot5 = plt.figure()
plt.plot(GS_serial_data['x'],GS_serial_data['y values'],'o',label='2(2-x^2-y^2)')
plt.plot(GS_serial_data['x'],GS_serial_data['y solved'],'x-',label='Gauss Seidel Serial')
plt.plot(GS_Diag_data_g01_thr2['x'],GS_Diag_data_g01_thr2['y solved'],'--',label='Gauss Seidel Diagonal Wave')
plt.plot(GS_BR_data_g01_thr2['x'],GS_BR_data_g01_thr2['y solved'],'-',label='Gauss Seidel Red Black')
plt.xlabel('x')
plt.ylabel('f(x,y)')
plt.legend()
plt.grid()
plt.title('Gauss Seidel comparison data for y=0.5')
fig5=plot5.get_figure()
fig5.savefig('Result_Images_and_Docs/Gauss_Seidel_comparison.png')

plot6 =plt.figure()
plt.plot(GS_DiagWave_grid_time['grid size'],GS_DiagWave_grid_time['time taken'],'-o' ,label='Diagonal Wave method')
plt.plot(GS_RedBlack_grid_time['grid size'],GS_RedBlack_grid_time['time taken'], '-x',label='Red Black method')
plt.xlabel('grid size')
plt.ylabel('time taken (in seconds)')
plt.title('Time taken for each grid size for number of threads = 8')
plt.legend()
plt.grid()
fig6=plot6.get_figure()
fig6.savefig('Result_Images_and_Docs/TimeVsGridSize_comparison.png')

plot7 =plt.figure()
plt.plot(GS_DiagWave_thread_time['grid size'],GS_DiagWave_thread_time['time taken'],'-o' ,label='Diagonal Wave method')
plt.plot(GS_RedBlack_thread_time['grid size'],GS_RedBlack_thread_time['time taken'], '-x',label='Red Black method')
plt.xlabel('grid size')
plt.ylabel('time taken (in seconds)')
plt.title('Time taken for number of threads for grid size = 0.005')
plt.legend()
plt.grid()
fig7=plot7.get_figure()
fig7.savefig('Result_Images_and_Docs/TimeVsNoOfThreads_comparison.png')

print('All images have been successfully saved in the workspace \n')