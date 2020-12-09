from __future__ import division
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np



def belonging_coeff(d, t, node):
    adj_mat = np.matrix([
                        [0, 1 / 3, 1 / 3, 1 / 3, 0, 0, 0, 0, 0, 0, 0],
                        [1 / 5, 0, 1 / 5, 1 / 5, 1 / 5, 1 / 5, 0, 0, 0, 0, 0],
                        [1 / 5, 1 / 5, 0, 1 / 5, 1 / 5, 0, 1 / 5, 0, 0, 0, 0],
                        [1 / 5, 1 / 5, 1 / 5, 0, 1 / 5, 0, 0, 1 / 5, 0, 0, 0],
                        [0, 1 / 9, 1 / 9, 1 / 9, 0, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9],
                        [0, 1 / 3, 0, 0, 1 / 3, 0, 1 / 3, 0, 0, 0, 0],
                        [0, 0, 1 / 4, 0, 1 / 4, 1 / 4, 0, 1 / 4, 0, 0, 0],
                        [0, 0, 0, 1 / 3, 1 / 3, 0, 1 / 3, 0, 0, 0, 0],
                        [0, 0, 0, 0, 1 / 3, 0, 0, 0, 0, 1 / 3, 1 / 3],
                        [0, 0, 0, 0, 1 / 3, 0, 0, 0, 1 / 3, 0, 1 / 3],
                        [0, 0, 0, 0, 1 / 3, 0, 0, 0, 1 / 3, 1 / 3, 0]
                        ])

    tran_mat = np.matrix([
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
                         ])

    def onestep(q):
        p = 0
        for i in G[q]:
            if i == 0:
                p = 1 / len(G[q])
        return p
    temp = tran_mat * adj_mat
    for k in G[node]:
        sum = (temp[k, 0]) / d

    sum = sum + onestep(node)
    for i in range(2, t):
        temp = temp * adj_mat
        for j in G[node]:
            temp2 = (temp[j, 0]) / d
            sum = sum + temp2
    return sum


def calculation_for_costFunction(params):

    n = len(params) - 1
    sum_i = n + (n*((n-1)/2))
    sum_i2 = 0
    sum_bc = 0
    sum_bc_i = 0

    for i in params:
        sum_i2 = sum_i2 + (i[0]**2)
        sum_bc = sum_bc + i[1]
        sum_bc_i = sum_bc_i + (i[0] * i[1])

    w = ((n*sum_bc_i) - (sum_i*sum_bc)) / ((n*sum_i2) - (sum_i**2))
    b = ((sum_i2*sum_bc) - (sum_i*sum_bc_i)) / ((n*sum_i2) - (sum_i**2))

    print (n, sum_i, sum_i2, sum_bc, sum_bc_i)
    return w, b

def cost_function(params, a, c, w, b):
    #print(params)

    p = min(a, c)
    q = max(a, c)
    cost_pq = 0
    for i in params[p:q+1]:
        eq = (w*i[0]) + b
        eq = eq - i[1]
        eq = eq**2
        cost_pq = cost_pq + eq
    return cost_pq

def rwlr(G):

    b_coeff = [(0, 0)]
    for i in range(1, 11):
        b_coeff.append((i, belonging_coeff(len(G[i]), 5, i)))
        print (b_coeff[i-1])

    def takefirst(ele):
        return ele[1]
    b_coeff.sort(key=takefirst, reverse=True)
    
    b_coeff = np.array(b_coeff)
    b_coeff_sort=b_coeff
    print(b_coeff_sort,"===:")

    F = b_coeff[:-1, 1]
    print(F,"=============================")
    Fu=b_coeff[:-1,0]
    print(Fu,"==========================")

    K = 1
    preQua = 0
    curQua = 0
    level = defaultdict(list)

    f = [[]]
    g = [[]]
    fi = []
    gi = []
    w, b = calculation_for_costFunction(b_coeff)
    print (w, b)
    for i in range(1, 11):
        cost1i = cost_function(b_coeff, 1, i, w, b)
        fi.append(cost1i)
        gi.append(-1)
    gi.insert(0, 0)
    fi.insert(0, 0)
    f.append(fi)
    g.append(gi)

    while True:

        preQua = curQua
        FL = []
        for fl in range(0, K):
            FL.append(F[fl])
        curQua = max(FL)

        if curQua < preQua or K == len(b_coeff) - 1: #or K == 3:
            break
        K = K + 1
        fi = []
        gi = []
        for i in range(1, 11):
            fki = []
            if i == 1:
                for j in range(1, K):
                    costki = f[K-1][j] + cost_function(b_coeff, j + 1, i, w, b)
                    fki.append(costki)
            if K == i:
                costki = f[K-1][K-1] + cost_function(b_coeff, K, i, w, b)
                fki.append(costki)
            else:
                for j in range(min(i-1, K-1), max(i, K)):
                    costki = f[K-1][j] + cost_function(b_coeff, j+1, i, w, b)
                    fki.append(costki)

            fi.append(min(fki))
            gi.append(np.argmin(fki) + 1)

        gi.insert(0, 0)
        fi.insert(0, 0)
        f.append(fi)
        g.append(gi)
    for x in f:
        print (x)
    for x in g:
        print (x)

        k = K
        i = len(b_coeff) - 1

        while True:
            if k == 0:
                break
            level[k] = []
            for j in range(min(i, g[k][i]), max(i, g[k][i])+1):
                level[k].append(j)
            i = g[k][i]
            k = k - 1
    K =  K - 1
    k = K
    i = len(b_coeff) - 1

    while True:
        if k == 0:
            break
        level[k] = []
        for j in range(min(i, g[k][i]), max(i, g[k][i])+1):
            level[k].append(j)
        i = g[k][i]
        k = k - 1
    for i in level:
        print (i, level[i])

    X = np.array(b_coeff)
    plt.scatter(X[:, 0], X[:, 1])
    plt.xlabel("Nodes")
    plt.ylabel("Coeffsss")
    plt.show()

adj_mat = [
    [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0],
    [0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0],
    [0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1],
    [0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0],
    [0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1],
    [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0],
]
G = defaultdict(list)
u, v = 0, 0
for i in adj_mat:
    u = 0
    for j in i:
        if j == 1:
            G[u].append(v)
        u += 1
    v += 1
print ("Graph...")
for i in range(len(G)+1):
    print (i, "->", G[i])
rwlr(G)