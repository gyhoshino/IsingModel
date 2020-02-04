import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

lower_limit = 2
upper_limit = 2.2
p0 = [100, 0, 3]

def f(x, a, b, critical_temp):
    return a*abs((x - critical_temp)/critical_temp)**(-b)
x = []
y = []
sigma = []
with open('grace_data/s0_N100/EM.csv','r') as csvfile:
#with open('data1/s1_N50/EM.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        try:
            x.append(abs(float(row[0])))
            y.append(abs(float(row[3])))
            sigma.append(float(row[6]))
        except:
            print("Error")
x = np.array(x)
y = np.array(y)
y = y**2/(x**2)
sigma = np.array(sigma)

index = []
for i in range(len(x)):
    if (x[i] > lower_limit) and (x[i] < upper_limit):
        index.append(i)
x = np.array(x[index])
y = np.array(y[index])
sigma = np.array(sigma[index])
full_fit = curve_fit(f, x, y, p0)
fit = full_fit[0]
var = full_fit[1]
FIT_a = fit[0]
FIT_b = fit[1]
FIT_CRITICAL_TEMP = fit[2]

error = ((2*y)/(x**2))*sigma

plt.errorbar(x,y, yerr=error, linewidth=1, marker='o', label="Data", linestyle='None')
plt.plot(x,f(x, FIT_a, FIT_b, FIT_CRITICAL_TEMP), label="Fit")
plt.xlabel('Temperature')
plt.ylabel('Heat Capacity')
plt.title('Heat Capacity vs Temperature')
plt.legend()
plt.show()

residuals = (y - f(x, FIT_a, FIT_b, FIT_CRITICAL_TEMP))**2
chi_squared = sum(residuals/error**2)/len(residuals)
print("alpha:" +  str(fit[1]))
print("Uncertainty alpha:" + str(np.sqrt(var[1,1])))
print("TC:" + str(fit[2]))
print("Uncertainty TC:" + str(np.sqrt(var[2,2])))
print("chi-squared:" + str(chi_squared))
