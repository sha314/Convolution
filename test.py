import convolution
import numpy as np
a = np.random.randint(0, 9, 9)
a = a.reshape((3,3))
a.dtype = np.float32
print(a)
b = convolution.convolute(a.tolist())
print(b)
