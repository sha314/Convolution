import convolution
import numpy as np
import matplotlib.pyplot as plt

filename = 'data_json.txt'
data_in = np.loadtxt(filename, usecols=(3,4))
# print(data_in)
data_out = convolution.convolute(data_in.tolist())
data_out = np.array(data_out)
p = np.loadtxt(filename, usecols=(0))
plt.plot(p, data_in[:,0], label='old')
plt.plot(p, data_out[:,0], label='new')
plt.legend()
plt.show()

plt.plot(p, data_in[:,1], label='old')
plt.plot(p, data_out[:,1], label='new')
plt.legend()
plt.show()