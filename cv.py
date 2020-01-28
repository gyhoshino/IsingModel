import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

lower_limit = 2
upper_limit = 2.2
p0 = [100, 0, 3]

def f(x, a, b, critical_temp):
    return a*abs((x - critical_temp)/critical_temp)**b
x = []
y = []
sigma = []
with open('data/s13_N100/EM.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        try:
            x.append(abs(float(row[0])))
            y.append(abs(float(row[2])))
            sigma.append(float(row[6]))
        except:
            print("Error")
x = np.array(x)
y = np.array(y)
sigma = np.array(sigma)

index = []
for i in range(len(x)):
    if (x[i] > lower_limit) and (x[i] < upper_limit):
        index.append(i)
x = np.array(x[index])
y = np.array(y[index])
sigma = np.array(sigma[index])
fit = curve_fit(f, x, y, p0)[0]
print(fit)
FIT_a = fit[0]
FIT_b = fit[1]
FIT_CRITICAL_TEMP = fit[2]

plt.plot(x,y, linewidth=1, marker='o', label="Data", linestyle='None')
plt.plot(x,f(x, FIT_a, FIT_b, FIT_CRITICAL_TEMP), label="Fit")
plt.xlabel('Temperature')
plt.ylabel('Heat Capacity')
plt.title('Heat Capacity vs Temperature')
plt.legend()
plt.show()

residuals = (y - f(x, FIT_a, FIT_b, FIT_CRITICAL_TEMP))**2
chi_squared = sum(residuals/sigma**2)/len(residuals)
print(chi_squared)