import csv
import matplotlib.pyplot as plt
import numpy as np

x = []
y = []

with open('data/s1_N20/EM.csv','r') as csvfile:
	plots = csv.reader(csvfile, delimiter=',')
	for row in plots:
		try:
			x.append(abs(float(row[0])))
			y.append(abs(float(row[4])))
		except:
			print("Error")

plt.plot(x,y, label='Magnetization', linewidth=1)
#plt.scatter(x,y, label='Magnetization')
plt.xlabel('Temperature')
plt.ylabel('Mean Magnetization')
plt.title('Magnetization vs Temperature')
plt.legend()
plt.show()