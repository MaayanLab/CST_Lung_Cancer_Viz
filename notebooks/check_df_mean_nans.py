# quick check of the behavior of ignore nans in dataframe mean

import numpy as np
import pandas as pd

zero_data = np.ones(shape=(10,2))

zero_data[:,1] = 2

# add some nans
zero_data[3,1] = np.NaN

zero_data[6,:] = np.NaN

zero_data[8,0] = np.NaN

tmp = pd.DataFrame(zero_data, columns=['A','B'])

print(tmp)

df_mean = tmp.mean(axis=1)
print('\n\n')
print(df_mean)