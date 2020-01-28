import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

E = []
for i in range(5, 205, 5):
    filename = 'data/s0_N' + str(i) + '/EM.csv'
    with open(filename) as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        for row in plots:
            try:
                E.append(abs(float(row[4])))
            except:
                print("Error")
    E = np.array(E)
print(E)