### Dataset is too small, so it is artifically increased ###

import pandas as pd
import numpy as np

mnt = "drive/MyDrive/HPC4DS/benchmarking/"

df = pd.read_csv(mnt+'Credit_Data_Light.csv', sep=',')
df.head()
print(df.shape)

# Fill missing values
df.isnull().any()
df['MINIMUM_PAYMENTS'] = df['MINIMUM_PAYMENTS'].fillna(df['MINIMUM_PAYMENTS'].median())
df.isnull().any()

# Specify the number of rows to append
n = 9000

# Create an empty DataFrame with the same columns as the original
new_rows = pd.DataFrame(columns=df.columns)

# Append n rows to the DataFrame
for i in range(n):
    new_row = {}
    for col in df.columns:
        min_val = df[col].min()
        max_val = df[col].max()
        new_row[col] = min_val + (max_val - min_val) * np.random.rand()
    new_rows = new_rows.append(new_row, ignore_index=True)

# Append the new rows to the original DataFrame
new_rows = df.append(new_rows, ignore_index=True)

print(new_rows.shape)
new_rows.tail()

new_rows.to_csv(mnt+'Credit_Data_Heavy.csv', index=False)
df.to_csv(mnt+'Credit_Data_Light.csv', index=False)