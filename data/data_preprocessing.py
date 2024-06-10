import pandas as pd
import numpy as np
import warnings
from sklearn.preprocessing import MaxAbsScaler, KBinsDiscretizer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split


def bin_data(data, max_bin):
    # Define the bin edges
    bins = np.linspace(0, 1, max_bin + 1)

    # Apply binning to each column
    binned_data = data.apply(lambda x: pd.cut(x, bins=bins, labels=False, include_lowest=True), axis=0)

    return binned_data




# Load the CSV file
file_path = './heart_failure_clinical_records.csv'
data = pd.read_csv(file_path)

# Drop the first five columns
data = data.iloc[:, 5:]

# Separate features and labels
features = data.iloc[:, :-1]  # All columns except the last one
labels = data.iloc[:, -1]  # The last column

# Normalize the features to range [0, 1] before binning
from sklearn.preprocessing import MinMaxScaler

scaler = MinMaxScaler()
normalized_features = scaler.fit_transform(features)
normalized_features_df = pd.DataFrame(normalized_features, columns=features.columns)

# Bin the normalized features
max_bin = 5  # Replace with your desired max_bin value
binned_features = bin_data(normalized_features_df, max_bin)

# Combine binned features and labels
binned_data = pd.concat([binned_features, labels], axis=1)

def custom_train_test_split(data, frac_of_training):
    # Check if 'target' column exists, if not assume it's the rightmost column
    if 'target' not in data.columns:
        target_column = data.columns[-1]
    else:
        target_column = 'target'

    # Separate the data based on the value of the "target" column
    target_0 = data[data[target_column] == 0]
    target_1 = data[data[target_column] == 1]

    # Split the group where the target is 0 into two parts: frac_of_training and remaining
    target_0_train, target_0_test = train_test_split(target_0, test_size=(1 - frac_of_training), random_state=42)

    # Combine the train part of the group where the target is 0 with all the rows where the target is 1
    train_set = target_0_train
    test_set = pd.concat([target_0_test, target_1])

    # Rename the target column
    train_set.rename(columns={target_column: 'target'}, inplace=True)
    test_set.rename(columns={target_column: 'target'}, inplace=True)

    # Save the resulting DataFrames into two separate CSV files with 'target' as the column name
    train_set.to_csv("train.csv", index=False)
    test_set.to_csv("test.csv", index=False)

custom_train_test_split(binned_data, 0.3)
