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
# csv_file = str(sys.argv[1])

# with open(csv_file,'r') as csvfile:
with open('Scaling.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        # if int(row[1])==32:
        x.append(int(row[0]))
        y.append(float(row[1]))
        z.append(float(row[2]))
        # x.append(math.log(int(row[0]),10))
        # y.append(math.log(float(row[1]),10))
        # z.append(math.log(float(row[2]),10))


# print(x, y, z)
plt.plot(x, y,label='Serial',color='blue')#,linestyle=':')
plt.plot(x, z,label='Parallel',color='red')#,linestyle=':')
# plt.plot(X1,Y1,label='Inc. M = 1024',color='blue')
# plt.plot(X1,Z1,label='Exc. M = 1024',color='red')
plt.xlabel('Iteration number (N)')
plt.ylabel('Program execution time (ms)')
plt.title('Program execution time vs. Iteration Number')
plt.subplots_adjust(bottom=.25, left=.25)
plt.legend(loc='best')
# plt.show()
plt.savefig('Scaling.png')
plt.close() 
