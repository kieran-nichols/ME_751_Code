import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
import sys
import numpy as np
import math

x = []
y = []
z = []
a = []
b = []
#c = []
# csv_file = str(sys.argv[1])

# with open(csv_file,'r') as csvfile:
i = 0
with open('coor_state.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        # if int(row[1])==32:
        x.append(i)
        y.append(float(row[0])) # horizontal position, x_pos
        z.append(float(row[3])) # forward vel, x_vel
        a.append(float(row[4]))	# vertical vel, y_vel
        #b.append(2)	# horizontal line as base ref for horz pos
        #a.append(int(row[0]))
        #b.append(float(row[1]))
        #c.append(float(row[2]))
        # x.append(math.log(int(row[0]),10))
        # y.append(math.log(float(row[1]),10))
        # z.append(math.log(float(row[2]),10))
        i+=1


# print(x, y, z)
plt.plot(x, y,label='x_pos (m)')
plt.plot(x, z,label='x_vel (m/s)')
m,c = np.polyfit(x,y,1)
x1 = np.array(x)
y_fit = m*x1+c
#print(y)
#print(type(np.array(x)),'\n',type(y_fit))
plt.plot(x1, y_fit,label='x_pos_fit')
#plt.plot(x, b,label='refernce line')
#plt.plot(x, a,label='y_vel')
plt.xlabel('Iteration number (N)')
plt.ylabel('Coor_state')
plt.title('Coordinate state vs. Iteration Number for serial code')
plt.subplots_adjust(bottom=.25, left=.25)
plt.legend(loc='best')
# plt.show()
plt.savefig('State.png')
plt.close() 
