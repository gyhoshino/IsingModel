import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy.special import ellipk

CRITICAL_TEMP = 2.269

T = []
average_E = []
average_M = []
std_E = []
std_M = []
stdstd_E = []
stdstd_M = []
cv = []
chi = []

with open('data/EM.csv','r') as csvfile:
	plots = csv.reader(csvfile, delimiter=',')
	for row in plots:
		try:
			if ((float(row[0]) > CRITICAL_TEMP) and (float(row[0]) < 2.5)):
				T.append(float(row[0]))
				average_E.append(float(row[2]))
				average_M.append(abs(float(row[4])))
				std_E.append(abs(float(row[3])))
				std_M.append(abs(float(row[5])))
				stdstd_E.append(abs(float(row[6])))
				stdstd_M.append(abs(float(row[7])))
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
stdstd_E = np.array(stdstd_E)
stdstd_M = np.array(stdstd_M)
cv = np.array(cv)
chi = np.array(chi)

def chi_squared(x, y, y_std, y_theor):
	 return sum((y - y_theor)**2 / (y_std)**2) / len(y)

def chi_fit(T, a, gamma, T_c):
	return a * abs((T - T_c) / T_c)**(-gamma)

def plot_chi_fit(log_plot = False):

	p0 = [1, 2, 2]
	a, gamma, T_c = curve_fit(chi_fit, T, chi, p0)[0]
	std_chi = (2 * std_M / T) * stdstd_M

	print("Gamma: ", gamma)
	print("Critical Temperature: ", T_c)
	print("Chi Squared: ", chi_squared(T, chi, std_chi, chi_fit(T, a, gamma, T_c)))

	if log_plot:
		plt.errorbar(np.log((T - T_c) / T_c), np.log(chi), yerr=(1 / chi) * std_chi, markersize=3, capsize=3, fmt='o', color='blue', zorder=1, label='Simulated data')
		plt.plot(np.log((T - T_c) / T_c), np.log(chi_fit(T, a, gamma, T_c)), color='red', zorder=0, label='Power law fit')
		plt.legend()
		plt.ylabel('Log Magnetic Susceptibility')
		plt.xlabel('Log Reduced Temperature')
		plt.title('Power Law Behavior of Magnetic Susceptibility | N = 100')
		plt.show()
	else:
		plt.errorbar(T, chi, yerr=std_chi, markersize=3, capsize=3, fmt='o', color='blue', zorder=1, label='Simulated data')
		plt.plot(T, chi_fit(T, a, gamma, T_c), color='red', zorder=0, label='Power law fit')
		plt.legend()
		plt.ylabel('Magnetic Susceptibility')
		plt.xlabel('Temperature (J/k)')
		plt.title('Power Law Behavior of Magnetic Susceptibility | N = 100')
		plt.show()


plot_chi_fit(log_plot=True)





