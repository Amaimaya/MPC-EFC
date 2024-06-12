import math
import time

import pandas as pd

# ------------ UTILS ------------

def unique_f(ar, return_index=False, return_inverse=False, return_counts=False):
    unique_list = []
    index_list = []
    inverse_list = [None] * len(ar)
    counts_dict = {}

    # Populate unique_list, index_list, and counts_dict
    for i, element in enumerate(ar):
        if element not in unique_list:
            unique_list.append(element)
            index_list.append(i)
            counts_dict[element] = 1
        else:
            counts_dict[element] += 1
        inverse_list[i] = unique_list.index(element)

    # Prepare the return_counts output
    if return_counts:
        counts_list = [counts_dict[element] for element in unique_list]

    # Prepare the output
    output = (unique_list,)
    if return_index:
        output += (index_list,)
    if return_inverse:
        output += (inverse_list,)
    if return_counts:
        output += (counts_list,)

    # Return output based on what the user has requested
    if len(output) == 1:
        return output[0]  # If only unique_list is requested, return it directly.
    return output


def equal_scalar(lst, scalar):
    return [element == scalar for element in lst]


def negative(values):
    if isinstance(values[0], list):  # Check if it's a multidimensional list
        return [[-x for x in row] for row in values]
    else:
        return [-x for x in values]


def log_2d(matrix):
    return [[math.log(element, 2.71828) for element in row] for row in matrix]


def log_1d(array):
    return [math.log(element, 2.71828) for element in array]


def cholesky_matrix_inverse(A):
    # https://github.com/marianpetruk/ad_fontes_2017/tree/master

    def cholesky(A):
        n = len(A)

        L = [[0.0] * n for _ in range(0, n, 1)]

        for i in range(0, n, 1):
            for j in range(0, i + 1, 1):
                tmp_sum = 0
                for k in range(0, j, 1):
                    tmp_sum += L[i][k] * L[j][k]
                if i == j:  # Diagonal elements
                    L[i][j] = math.sqrt(A[i][i] - tmp_sum)
                else:
                    L[i][j] = (1.0 / L[j][j] * (A[i][j] - tmp_sum))
        return L

    def transpose(a):
        n = len(a)
        matrix_transpose = a
        for i in range(0, n, 1):
            for j in range(0, i + 1, 1):
                matrix_transpose[j][i] = a[i][j]
        return matrix_transpose

    def intermediate_diagonal_inverse_matrix(r):
        n = len(r)
        S = [[0] * n for i in range(1, n + 1, 1)]
        for i in range(0, n, 1):
            S[i][i] = 1 / r[i][i]
        return S

    def krishnamoorthy_menon_method(R, S):
        n = len(R)
        X = [[0.0] * n for i in range(0, n, 1)]

        for i in range(n - 1, -1, -1):
            for j in range(i, -1, -1):
                tmp_sum = 0
                for k in range(j + 1, n, 1):
                    tmp_sum += R[j][k] * X[k][i]
                X[j][i] = (S[j][i] - tmp_sum) / R[j][j]
                X[i][j] = X[j][i]
        return X

    L = cholesky(A)
    R = transpose(L)
    S = intermediate_diagonal_inverse_matrix(R)
    X = krishnamoorthy_menon_method(R, S)
    return X

def zeros_1d(shape_a):
    return [0 for _ in range(shape_a)]

def zeros_2d(shape_a, shape_b):
    return [[0 for _ in range(shape_b)] for _ in range(shape_a)]

def zeros_3d(shape_a, shape_b, shape_c):
    return [[[0 for _ in range(shape_c)] for _ in range(shape_b)] for _ in range(shape_a)]

def zeros_4d(shape_a, shape_b, shape_c, shape_d):
    return [[[[
        0 for _ in range(shape_d)]
        for _ in range(shape_c)]
        for _ in range(shape_b)]
        for _ in range(shape_a)]

def dot_product(x, y):
    dp = 0
    for i in range(len(x)):
        dp += (x[i]*y[i])
    return dp

def dot_product_matrix(x, y):
    dp = 0
    for i in range(len(x)):
        dp += dot_product(x[i], y[i])
    return dp

def exponential_matrix(input_matrix):
    return [[2.71828 ** element for element in row] for row in input_matrix]


def multiply_2d_matrix_by_scalar(matrix, scalar):
    return [[element * scalar for element in row] for row in matrix]

def add_2d_matrix_by_scalar(matrix, scalar):
    return [[element + scalar for element in row] for row in matrix]

def multiply_4d_matrix_by_scalar(matrix, scalar):
    return [[[[element * scalar for element in d] for d in c] for c in b] for b in matrix]

def add_4d_matrix_by_scalar(matrix, scalar):
    return [[[[element + scalar for element in d] for d in c] for c in b] for b in matrix]


def divide_2d_matrix_by_scalar(matrix, scalar):
    return [[element / scalar for element in row] for row in matrix]


def divide_4d_matrix_by_scalar(matrix, scalar):
    return [[[[element / scalar for element in d] for d in c] for c in b] for b in matrix]

def cantor(x, y):
    output = [0.0] * len(x)
    for i in range(len(x)):
        output[i] = (x[i] + y[i]) * (x[i] + y[i] + 1) / 2 + y[i]
    return output

# ------------------ EFC ------------------

def site_freq(X_view, psdcounts, max_bin):
    print("!!!!!!!!!!!!!!!")
    print(X_view, psdcounts, max_bin)
    print("!!!!!!!!!!!!!!!")
    n_attr = len(X_view[1])  # total feature
    sitefreq = zeros_2d(n_attr, max_bin)

    for i in range(n_attr):
        for aa in range(max_bin):
            sitefreq[i][aa] = sum(equal_scalar([row[i] for row in X_view], aa))
    
    print("Site freq 1:", sitefreq)

    sitefreq = divide_2d_matrix_by_scalar(sitefreq, len(X_view))

    sitefreq = multiply_2d_matrix_by_scalar(sitefreq, 1 - psdcounts)

    sitefreq = add_2d_matrix_by_scalar(sitefreq, (psdcounts / max_bin))

    print("Site freq: ", sitefreq)
    return sitefreq

def pair_freq(X_view, sitefreq_view, psdcounts, max_bin):
    n_inst, n_attr = len(X_view), len(X_view[1])
    pairfreq = zeros_4d(n_attr, max_bin, n_attr, max_bin)

    for i in range(n_attr):
        for j in range(n_attr):
            c = cantor([row[i] for row in X_view], [row[j] for row in X_view])
            unique, aaIdx = unique_f(c, return_index=True)
            for x, item in enumerate(unique):
                pairfreq[i][X_view[aaIdx[x]][i]][j][X_view[aaIdx[x]][j]] = sum(equal_scalar(c, item))

    pairfreq = divide_4d_matrix_by_scalar(pairfreq, n_inst)
    pairfreq = multiply_4d_matrix_by_scalar(pairfreq, 1 - psdcounts)
    pairfreq = add_4d_matrix_by_scalar(pairfreq, psdcounts / (max_bin ** 2))

    for i in range(n_attr):
        for ai in range(max_bin):
            for aj in range(max_bin):
                if ai == aj:
                    pairfreq[i][ai][i][aj] = sitefreq_view[i][ai]
                else:
                    pairfreq[i][ai][i][aj] = 0.0

    print("PAIRFREQ: ", pairfreq)
    return pairfreq

def coupling(pairfreq_view, sitefreq_view, max_bin):
    n_attr = len(sitefreq_view)
    corr_matrix = zeros_2d(n_attr * (max_bin - 1), n_attr * (max_bin - 1))


    for i in range(n_attr):
        for j in range(n_attr):
            for ai in range(max_bin - 1):
                for aj in range(max_bin - 1):
                    corr_matrix[i * (max_bin - 1) + ai][j * (max_bin - 1) + aj] = (
                            pairfreq_view[i][ai][j][aj] - sitefreq_view[i][ai] * sitefreq_view[j][aj])

    print("Corr matrix:", corr_matrix)
    inv_corr = cholesky_matrix_inverse(corr_matrix)
    negative_m = negative(inv_corr)
    result = exponential_matrix(negative_m)
    return result

def local_fields(coupling_view, sitefreq_view, max_bin):
    n_inst = len(sitefreq_view)
    fields = [0.0] * (n_inst * (max_bin - 1))

    for i in range(n_inst):
        for ai in range(max_bin - 1):
            fields[i * (max_bin - 1) + ai] = sitefreq_view[i][ai] / sitefreq_view[i][max_bin - 1]
            acc = 1
            for j in range(n_inst):
                for aj in range(max_bin - 1):
                    acc *= (
                            coupling_view[i * (max_bin - 1) + ai][j * (max_bin - 1) + aj] ** sitefreq_view[j][aj])
            
            print("acc: ", acc)
            fields[i * (max_bin - 1) + ai] /= acc
    print("Localfields", fields)
    return fields

def compute_energy(X, coupling_matrix, local_fields, max_bin):
    n_inst, n_attr = len(X), len(X[0])
    energies = zeros_1d(n_inst)

    for i in range(n_inst):
        e = 0
        for j in range(n_attr):
            j_value = X[i][j]  # Get the value of attribute j for the current instance i
            if j_value != (max_bin - 1):
                for k in range(j, n_attr):
                    k_value = X[i][k]
                    if k_value != (max_bin - 1): # Check if the attribute value is not the last bin
                        # The "−∑e" part of the equation
                        e -= (coupling_matrix[j * (max_bin - 1)
                                                    + j_value][k * (max_bin - 1) + k_value])
                # The "−∑h" part of the equation
                e -= (local_fields[j * (max_bin - 1) + j_value])
        energies[i] = e
    
    print("energies", energies)
    return energies


def fit(X, pseudocounts, max_bin):
    sitefreq = site_freq(X, pseudocounts, max_bin)
    pairfreq = pair_freq(X,sitefreq, pseudocounts, max_bin)
    coupling_matrix = coupling(pairfreq, sitefreq, max_bin)
    localfields = local_fields(coupling_matrix, sitefreq, max_bin)
    coupling_matrix = log_2d(coupling_matrix)
    localfields = log_1d(localfields)

    print("localfields base", localfields)
    print("coupling matrix base", coupling_matrix)

    cutoff = define_cutoff(X, 0.95, coupling_matrix, localfields, max_bin)
    return sitefreq, pairfreq, coupling_matrix, localfields, cutoff

def define_cutoff(X, cutoff_quantile, coupling_matrix, local_fields, max_bin):
    energies = compute_energy(X, coupling_matrix, local_fields, max_bin)
    energies = sorted(energies)
    print("cutoff", energies[int(len(energies) * cutoff_quantile)])
    return energies[int(len(energies) * cutoff_quantile)]

def predict(X, sitefreq, pairfreq, coupling_matrix, localfields, cutoff, max_bin):
    energies = [compute_energy(X, coupling_matrix, localfields, max_bin)]
    y_pred = []
    for row in range(len(X)):
        minimum_energy = min([energies[i][row] for i in range(len(energies))])
        if minimum_energy < cutoff:
            y_pred.append(0)
        else:
            y_pred.append(1)
    return y_pred

def calculate_metrics(y_test, y_pred, positive_label=1):
    true_positive = 0
    false_positive = 0
    false_negative = 0

    for y_t, y_p in zip(y_test, y_pred):
        if y_p == positive_label and y_t == positive_label:
            true_positive += 1
        elif y_p == positive_label and y_t != positive_label:
            false_positive += 1
        elif y_p != positive_label and y_t == positive_label:
            false_negative += 1

    precision = true_positive / (true_positive + false_positive) if true_positive + false_positive != 0 else 0.0
    recall = true_positive / (true_positive + false_negative) if true_positive + false_negative != 0 else 0.0

    return precision, recall

 # ------------------ TESTING ------------------

X_train = pd.read_csv("../data/mediumbig/train.csv")
X_test = pd.read_csv("../data/mediumbig/test.csv")

y_train = X_train["target"].values.tolist()
y_test = X_test["target"].values.tolist()

X_train = X_train.drop(columns=["target"]).values.tolist()
X_test = X_test.drop(columns=["target"]).values.tolist()

max_bin = max(max(max(row) for row in X_train), max(max(row) for row in X_test)) + 1

start = time.time()
sitefreq, pairfreq, coupling_matrix, localfields, cutoff = fit(X_train, 0.5, max_bin)
stop = time.time()
print("Fit: ", stop - start, " seconds")
start = time.time()
y_pred = predict(X_test, sitefreq, pairfreq, coupling_matrix, localfields, cutoff,  max_bin)
stop = time.time()
print("Predict: ", stop - start, " seconds")

#print("y pred", y_pred)

precision, recall = calculate_metrics(y_test, y_pred)

#print("Precision:", precision)
#print("Recall:", recall)
