import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.special import ellipk

CRITICAL_TEMP = 2.269

def chi_squared(x, y, y_std, y_theor):
	 return sum((y-y_theor)**2 / (y_std)**2) / len(y)


T = []
average_E = []
average_M = []
std_E = []
std_M = []
cv = []
chi = []

with open('data/s4_N50/EM.csv','r') as csvfile:
	plots = csv.reader(csvfile, delimiter=',')
	for row in plots:
		try:
			if (float(row[0]) != 0):
				if (float(row[3]) != 0):
					T.append(float(row[0]))
					average_E.append(float(row[2]))
					average_M.append(abs(float(row[4])))
					std_E.append(abs(float(row[3])))
					std_M.append(abs(float(row[5])))
					# Specific heat in units of 1/k_b
					cv.append(float(row[3])**2 / float(row[0])**2)
					# Magnetic susceptibility in units of 1 / k_b
					chi.append(float(row[5])**2 / float(row[0]))
		except:
			pass


# Convert to numpy arrays
T = np.array(T)
E = np.array(average_E)
M = np.array(average_M)
std_E = np.array(std_E)
std_M = np.array(std_M)
cv = np.array(cv)
chi = np.array(chi)

#x = np.arange(0.01, 5, 0.0001)
x = T
analytic_E = -2*np.tanh(2/x) - (np.sinh(2/x)**2 - 1) / (np.sinh(2/x) * np.cosh(2/x)) * (2 * ellipk((2 * np.sinh(2/x) / np.cosh(2/x)**2)**2) / np.pi - 1)

analytic_M_1 = np.array([(1 - np.sinh(2/t)**(-4))**(1/8) for t in T if t < CRITICAL_TEMP])
analytic_M_2 = np.array([0 for t in T if t > CRITICAL_TEMP])
analytic_M = np.append(analytic_M_1, analytic_M_2)



"""
x2 = np.arange(0.01, CRITICAL_TEMP, 0.001)
analytic_M = (1 - np.sinh(2/x2)**(-4))**(1/8)
x2 = np.append(x2, CRITICAL_TEMP)
analytic_M = np.append(analytic_M, 0)
x3 = np.arange(CRITICAL_TEMP, 5, 0.01)
zero = 0*x3


plt.errorbar(T, M, yerr=std_M, markersize=3, capsize=3, fmt='o', color='blue', zorder=0, label='Simulated data')
plt.plot(x2, analytic_M, color='red', zorder=1, label='Theoretical curve')
plt.plot(x3, zero, color='red', zorder=2)
#plt.plot(x, analytic_E, color='red', zorder=1, label='Theoretical curve')
plt.xlabel('Temperature (J/k)')
plt.ylabel('Mean Magnetization')
plt.title('Magnetization vs. Temperature | N = 50')
plt.legend()
plt.show()
"""

print(chi_squared(T, M, std_M, analytic_M))



"""
# Reduced temperature calculation
reduced_t = np.array([abs((t - CRITICAL_TEMP) / CRITICAL_TEMP) for t in T if t < CRITICAL_TEMP])
log_reduced_t = np.log(reduced_t)

plt.errorbar(np.log(reduced_t), np.log(chi[:len(reduced_t)]), markersize=3, capsize=3, fmt='o', color='blue', zorder=0)
#plt.plot(x, intercept + slope * x, color='red', zorder=1)
plt.xlabel('Temperature')
plt.ylabel('Magnetic Susceptibility')
plt.title('N = 40  |  R Flip: 0.04  |  Annealing Steps: 80k | Burning Steps: 50k')
plt.show()


# Power law for SPECIFIC HEAT
reduced_t = np.array([abs((t - CRITICAL_TEMP) / CRITICAL_TEMP) for t in T if t < CRITICAL_TEMP])
log_reduced_t = np.log(reduced_t)
log_cv = np.log(cv[:len(reduced_t)])

# Regression on specific heat
slope, intercept, r_value, p_value, std_err = stats.linregress(log_reduced_t, log_cv)

print("Slope: %f" % slope)
print("r-value: %f" % p_value)


x = np.arange(-8, -3, 0.01)


#plt.errorbar(T, average_M, yerr=std_M, label='Magnetization', markersize=3, capsize=3, fmt='o', color='blue', zorder=0)
plt.errorbar(log_reduced_t, log_cv , markersize=3, capsize=3, fmt='-o', color='blue', zorder=0)
plt.plot(x, intercept + slope * x, color='red', zorder=1)
plt.xlabel('Log of Reduced Temperature')
plt.ylabel('Log of Specific Heat')
plt.title('N = 40  |  R Flip: 0.04  |  Annealing Steps: 80k | Burning Steps: 50k')
plt.show()

"""
