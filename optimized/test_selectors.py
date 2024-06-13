import time
import pandas as pd
import numpy as np

SIZE = "small"

print("Running test selectors...")

X = pd.read_csv("../data/" + SIZE + "/test.csv")
X_train = pd.read_csv("../data/" + SIZE + "/train.csv")
max_value_X = X.max().max()
max_value_X_train = X_train.max().max()
max_bin = max(max_value_X, max_value_X_train) + 1
n_inst = X_train.shape[0]
n_attr = X_train.shape[1]-1 # because +1 due to y_pred
len_local_fields = n_attr * (max_bin - 1)
len_coupling_matrix = n_attr * (max_bin - 1)
def compute_selectors():
    all_fields_selectors = []
    all_matrix_selectors = []

    def zeros_1d(shape_a):
        return [0 for _ in range(shape_a)]

    def zeros_2d(shape_a, shape_b):
        return [[0 for _ in range(shape_b)] for _ in range(shape_a)]
    
    n_inst = X.shape[0] 
    n_attr = X.shape[1] - 1

    for i in range(n_inst):
        matrix_selector1 = zeros_1d(len_local_fields)
        matrix_selector2 = zeros_2d(len_coupling_matrix, len_coupling_matrix)

        for j in range(n_attr):
            j_value = X.iloc[i, j]
            is_j_value = 1 if (j_value != (max_bin-1)) else 0

            for k in range(j, n_attr):
                k_value = max_bin-1 if is_j_value == 0 else X.iloc[i, k]
                is_k_value = 1 if (k_value != (max_bin-1)) else 0

                stop = False
                for l in range(len_coupling_matrix):
                    if stop:
                        break
                    for m in range(len_coupling_matrix):
                        if l == (j * (max_bin - 1) + j_value) and m == (k * (max_bin - 1) + k_value) and is_k_value:
                            matrix_selector2[l][m] = 1
                            stop = True
                            break

            for k in range(len_local_fields):
                if k == (j * (max_bin - 1) + j_value) and is_j_value:
                    matrix_selector1[k] = 1
                    break

        all_fields_selectors.append(matrix_selector1)
        all_matrix_selectors.append(matrix_selector2)
    
    return np.array(all_fields_selectors), np.array(all_matrix_selectors)

# Compute selectors
start = time.time()
all_fields_selectors, all_matrix_selectors = compute_selectors()
stop = time.time()
print("Test selector took: ", stop - start, " seconds")

# Save arrays to numpy files
np.save("../data/" + SIZE + "/fields_selectors.npy", all_fields_selectors)
np.save("../data/" + SIZE + "/matrix_selectors.npy", all_matrix_selectors)

print("Finished running test selectors...")