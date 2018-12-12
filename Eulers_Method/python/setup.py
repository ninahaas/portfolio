## Nina Haas
## This python program calls the makefile to compile the fortran code
## Then plots the output and calculates the error

import numpy as np
import os
import matplotlib.pyplot as plt

error = np.zeros((4))


## real_sol()
##Input: time
##Output: real solution to ODE
def real_sol(t):
	return -(2*np.log(t**2 + 1) + 4)**(0.5)

## create_data()
## Input: output text file, number of lines
## Output: matrix with rows being each data points
def create_data(filename,num):
	data = np.zeros((num,4))
	f = open(filename,'r')
	i = 0
	for line in f:
        	data[i,0] = float(line.split()[0])      # time value
        	data[i,1] = float(line.split()[1])      # solution
		data[i,2] = real_sol(data[i][0])        # real solution
		data[i,3] = abs(data[i,2] - data[i,1])  # error
        	i = i+1
	f.close()
	return data


## create_plot()
## Input: data matrix, total error, output figure name
## Output: None
def create_plot(data, error, outname):
	x = np.arange(len(data))
	plt.figure()
	plt.plot(x, data[:,1], color='red', marker='o', linestyle='dashed', linewidth=1, markersize=4)
	plt.xticks(fontsize = 7)
	plt.yticks(fontsize = 7)
	plt.plot(x, data[:,2], color='blue', linestyle='solid', linewidth=1)
	plt.grid(color='gray',linestyle='-')
	plt.xlabel('t axis', fontsize=10)
	plt.ylabel('y axis', fontsize=10)
	plt.title(error)
	plt.savefig(outname)


## MAIN ############################################################


## run makefile to compile fortran code
os.chdir('../fortran')
os.system('make clean')
os.system('make')
os.system('./main.exe')

## create data matrix for each output
data_8 = create_data('output_8.txt',8)
data_16 = create_data('output_16.txt',16)
data_32 = create_data('output_32.txt',32)
data_64 = create_data('output_64.txt',64)

## calculate total error of each output
error[0] = np.sum(data_8[:,3])
error[1] = np.sum(data_16[:,3])
error[2] = np.sum(data_32[:,3])
error[3] = np.sum(data_64[:,3])

os.chdir('../python')

## create and save plots to python folder
create_plot(data_8, error[0], 'result_8.png')
create_plot(data_16, error[1], 'result_16.png')
create_plot(data_32, error[2], 'result_32.png')
create_plot(data_64, error[3], 'result_64.png')
