#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from matplotlib import dates as mpl_dates

# Read data
data = np.loadtxt("history.txt", dtype=int)

# Split columns
date = []
for i in range(0, data.shape[0]):
   date.append(datetime(data[i,0], data[i,1], data[i,2], data[i,3], data[i,4], data[i,5]))
fortran = data[:,6]
cpp = data[:,7]
python = data[:,8]

# Minimum date
imin = 0
for i in range(0, data.shape[0]):
   if fortran[i]+cpp[i]+python[i] > 0:
      imin = i

# Plot data
fig, ax = plt.subplots()
ax.plot(date,fortran,'-r')
ax.plot(date,cpp,'-b')
ax.plot(date,python,'-g')

# X axis format
plt.gcf().autofmt_xdate()
date_format = mpl_dates.DateFormatter('%d-%m-%Y')
plt.gca().xaxis.set_major_formatter(date_format)
plt.xlim(xmin=date[imin])
ax.set_xlabel("Date")

# Y axis format
plt.ylim(ymin=0,ymax=1.1*max(fortran))
ax.set_ylabel("Number of lines")

# Title
plt.title('SABER evolution')
plt.ylabel('Number of lines')
plt.tight_layout()

# Save figure
plt.savefig("history.jpg", format="jpg", dpi=300)
plt.close()
