# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 12:32:28 2021

@author: jdiab
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

t1=[-20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
s7p=[14.1179, 13.2404, 12.3614, 11.5127, 10.7202, 10.0007, 9.3615, 8.8020, 8.3170, 7.8989, 7.5390, 7.2293] 
s5p=[4.8141, 4.4541, 4.1635, 3.9260, 3.7300, 3.5666, 3.4293, 3.3131, 3.2142, 3.1294, 3.0565, 2.9934]
sns.set()

plt.figure(1)
plt.plot(t1, s7p, label='pH 7')
plt.plot(t1, s5p, label='pH 5 (unstable)')
plt.xlabel("Temperature (\N{DEGREE SIGN}C)")
plt.ylabel("Size (nm)")
plt.title("Temp vs Size of Silica Colloids in Pure Water")
plt.legend()