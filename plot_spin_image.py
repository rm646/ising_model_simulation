## short script to plot spins image
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv(
    'spin_image_data.tsv',
    sep='\t',
    header=None,
    names=['x', 'y', 'b'],
    skiprows=6
)

# Get matrix size
x_range = [df.x.min(), df.x.max()]
y_range = [df.y.min(), df.y.max()]

if x_range[0]!=0 or y_range[0]!=0:
    raise ValueError('Bad image data: missing 0 index for x or y.')

spin_array = np.zeros((x_range[1]+1, y_range[1]+1))

df_ones = df[df.b == 1]
for index, row in df_ones.iterrows():
    spin_array[row.x, row.y] = row.b

plt.matshow(spin_array)
plt.show()
