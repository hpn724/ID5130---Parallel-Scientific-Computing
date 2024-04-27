import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np

def plot_3d_surface(dataframe):
    # Extract x and y values from the DataFrame
    x_values = dataframe.columns.astype(float)
    y_values = dataframe.index.astype(float)
    
    # Create a meshgrid for x and y values
    x_mesh, y_mesh = np.meshgrid(x_values, y_values)
    
    # Convert DataFrame values to a 2D array
    z_values = dataframe.values
    
    # Create a 3D figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the 3D surface
    ax.plot_surface(x_mesh, y_mesh, z_values, cmap='viridis')
    
    # Set labels and title
    ax.set_xlabel('Distance')
    ax.set_ylabel('Time')
    ax.set_zlabel('Value')
    ax.set_title('3D Surface Plot')
    
    # Show the plot
    plt.show()


Poisson_Jacobi_MPI_2_0_1 = pd.read_csv('Solution/Q1b/1D_Wave_MPI_2_upwind.csv',index_col=0,header=0)

plot_3d_surface(Poisson_Jacobi_MPI_2_0_1)