import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy.special import ellipk

CRITICAL_TEMP = 2.269
MAX_DISTANCE = 40 # 40 ORIGINALLY

T_array = []
spin_C_array = []
std_spin_C_array = []
xi_array = []
std_xi_array = []

with open('data/SpinCorr.csv','r') as csvfile:
	plots = csv.reader(csvfile, delimiter=',')
	for row in plots:
		try:
			if ((float(row[0]) >= 2) and (float(row[0]) <= 2.5)):
				T_array.append(float(row[0]))
				spin_C_array.append(np.array(list(float(row[i]) for i in range(2, 2 + MAX_DISTANCE))))
				std_spin_C_array.append(np.array(list(float(row[i]) for i in range(2 + 49, 2 + 49 + MAX_DISTANCE))))
		except:
			pass

distance = np.array([i for i in range(1, MAX_DISTANCE + 1)])
spin_C = np.array(spin_C_array)
std_spin_C = np.array(std_spin_C_array)

def chi_squared(x, y, y_std, y_theor):
	 return sum((y - y_theor)**2 / (y_std)**2) / len(y)

def C(x, xi, a, b):
	return a * np.exp(-x / xi) + b

def fit_C(distance, T, spin_C, std_spin_C, plot = False):

	p0 = [1, 2, 0]
	params = curve_fit(C, distance, spin_C, p0)
	xi, a, b = params[0]
	covariance = params[1]
	std_xi = np.sqrt(covariance[0,0])

	print("xi: ", xi)
	#print("a: ", a)
	print("std_xi: ", std_xi)

	chi_error = chi_squared(distance, spin_C, std_spin_C, C(distance, xi, a, b))
	print("Chi squared: ", chi_error)


	x = np.arange(1, 40, 0.01)

	if plot:
		plt.errorbar(distance, spin_C, yerr=std_spin_C, markersize=3, capsize=3, fmt='o', color='blue', zorder=1, label='Simulated data')
		plt.plot(x, C(x, xi, a, b), color='red', zorder=0, label='Fit')
		plt.legend()
		plt.ylabel('Correlation Function (T = %.2f)' % T)
		plt.xlabel('Distance (Adjacent Spins)')
		plt.title('Spin Correlation vs. Distance | N = 100')
		plt.show()

	return xi, std_xi

fit_C(distance, T_array[50], spin_C_array[50], std_spin_C_array[50], plot = True)

def xi_func(T, a, nu, T_c):
	return a * abs((T - T_c) / T_c)**(-nu)


def fit_xi(plot = True):
	
	for i in range(0, len(T_array)):
		xi, std_xi = fit_C(distance, T_array[i], spin_C_array[i], std_spin_C_array[i])
		xi_array.append(xi)
		std_xi_array.append(std_xi)

	T_above = []
	xi_above = []
	std_xi_above = []
	for i in range(len(T_array)):
		if (T_array[i] >= CRITICAL_TEMP):
			T_above.append(T_array[i])
			xi_above.append(xi_array[i])
			std_xi_above.append(std_xi_array[i])

	p0 = [1, 2, 2]
	params = curve_fit(xi_func, T_above, xi_above, p0)
	a, nu, T_c = params[0]
	covariance = params[1]
	std_nu = np.sqrt(covariance[1, 1])
	std_T_c = np.sqrt(covariance[2, 2])
	print("Nu: ", nu)
	print("T_c: ", T_c)
	#print("a: ", a)
	print("std Nu: ", std_nu)
	print("std T_c: ", std_T_c)

	x = np.arange(T_above[0], T_above[-1], 0.001)

	if plot:
		plt.errorbar(T_above, xi_above, yerr=std_xi_above, markersize=3, capsize=3, fmt='o', color='blue', zorder=1, label='Simulated data')
		plt.plot(x, xi_func(x, a, nu, T_c), color='red', zorder=0, label='Fit')
		plt.legend()
		plt.ylabel('Xi (Correlation Length)')
		plt.xlabel('Temperature (J/k)')
		plt.title('Correlation Length vs. Temperature | N = 100')
		plt.show()

fit_xi()










