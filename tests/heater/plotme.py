import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('./work/plant_res.csv', sep=",")

# Create the figure and axes
fig, ax = plt.subplots(nrows=1,ncols=1)
ax.plot(df["time"], df["htr.Tf"], label='fluid', color='blue', linestyle='-')
ax.plot(df["time"], df["htr.Tw"], label='wall', color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Temp. (C)')
ax.set_title('Heater/Cooler Test')
ax.legend()
ax.grid(True) # Add a grid for better readability
plt.show()
