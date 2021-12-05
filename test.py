import numpy as np
import matplotlib.pyplot as plt
from pylab import *


x = np.random.random(10)
a = 1.0
u = x * a / (1.+a)-1
cos_theta = (2. / u + 2.)/a + 1
print(cos_theta)
sin_theta = np.sqrt(1. - cos_theta**2)
print(sin_theta)