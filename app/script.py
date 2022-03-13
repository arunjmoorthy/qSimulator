import numpy
import math
import scipy
import scipy.linalg
import random

def trdist(A, B):
    U, D, V = numpy.linalg.svd(numpy.subtract(A, B))
    return numpy.sum(D)

iMat = [[1, 0], [0, 1]]
xMat = [[0, 1], [1, 0]]
yMat = [[0, complex(0, -1)], [complex(0, 1), 0]]
zMat = [[1, 0], [0, -1]]

def matChoose(ch):
    if ch == 'i':
        return iMat
    elif ch == 'x':
        return xMat
    elif ch == 'y':
        return yMat
    else:
        return zMat

def kronFor(paulis):
    if len(paulis) == 1:
        return numpy.matrix(matChoose(paulis[0]))
    paulis = str(paulis[::-1])
    res = numpy.kron(matChoose(paulis[1]), matChoose(paulis[0]))

    for i in range(2, len(paulis)):
        res = numpy.kron(matChoose(paulis[i]), res)
    return res


def TrotMoreArgs(ops, x):
    res = numpy.identity(len(ops[0]))
    for i in range(len(ops)):
        mat = ops[i]
        mat = mat * x * 0.5
        fin = scipy.linalg.expm(mat)
        res = res.dot(fin)
    for i in range(len(ops) - 1, -1, -1):
        mat = ops[i]
        mat = mat * x * 0.5
        fin = scipy.linalg.expm(mat)
        res = res.dot(fin)
    return res

def idealexp(ops, x):
    res = numpy.matrix(ops[0])
    for i in range(1, len(ops)):
        res = numpy.add(res, numpy.matrix(ops[i]))
    return scipy.linalg.expm(x * res)

def recurse(m2, x, ops):
    if (m2 <= 2):
        return TrotMoreArgs(ops, x)
    p = pow((4 - pow(4, (1. / (m2 - 1)))), -1)
    fin = (numpy.linalg.matrix_power(recurse((m2 - 2), (p * x), ops), 2).dot(
        recurse((m2 - 2), ((1 - (4 * p)) * x), ops))).dot(numpy.linalg.matrix_power(recurse((m2 - 2), (p * x), ops), 2))
    return fin

def comm(p1, p2):
    comm = 0
    n = int(len(p1) / 2)
    for i in range(n):
        comm = comm + (p1[i] * p2[i + n] - p1[i + n] * p2[i])
    return comm % 2

def sympmatrix(paulis):
    m = numpy.zeros([len(paulis), len(paulis)])
    for i in range(len(m)):
        for j in range(len(m)):
            m[i, j] = comm(paulis[i], paulis[j])
    return m

def symptoadj(symp):
    for i in range(len(symp)):
        for j in range(len(symp)):
            if (not (i == j)):
                symp[i, j] = 1 - symp[i, j]
    return symp.tolist()

def find_clique(adj):
    clique = []
    rvertex = random.randrange(0, len(adj), 1)
    clique.append(rvertex)
    for v in range(len(adj)):
        if v in clique:
            continue
        toadd = True
        for u in clique:
            if adj[v][u] == 1:
                continue
            else:
                toadd = False
                break
        if (toadd):
            clique.append(v)
    return sorted(clique)

def remove_verts(adj, vrem):
    vrem = sorted(vrem, reverse=True)
    for v in vrem:
        for r in adj:
            del r[v]
    for v in vrem:
        del adj[v]
    return adj

def recover_labels(parts, totalv):
    verts = []
    for i in range(totalv):
        verts.append(i)
    labels = []
    for part in parts:
        part = sorted(part, reverse=True)
        piece = []
        for index in part:
            piece.append(verts[index])
        for index in part:
            del verts[index]
        labels.append(piece)
    return labels

def partition(adj):
    rolling = adj
    parts = []
    while (len(rolling) > 0):
        clique = find_clique(rolling)
        parts.append(clique)
        remove_verts(rolling, clique)
    return parts

def randompaulis(v, n):
    paulis = []
    while (len(paulis) < v):
        pauli = []
        for i in range(2 * n):
            pauli.append(random.randrange(0, 2))
        if (not (pauli in paulis)):
            paulis.append(pauli)
    return paulis

def symptop(psymp):
    pauli = ""
    n = int(len(psymp) / 2)
    for i in range(n):
        if (psymp[i] == 0 and psymp[i + n] == 0):
            pauli = pauli + "i"
        elif (psymp[i] == 1 and psymp[i + n] == 0):
            pauli = pauli + "x"
        elif (psymp[i] == 0 and psymp[i + n] == 1):
            pauli = pauli + "z"
        else:
            pauli = pauli + "y"
    return pauli

def clique_matrices(rpaulis):
    out = []
    adjmat = symptoadj(sympmatrix(rpaulis))
    parts = recover_labels(partition(adjmat), len(rpaulis))
    for i in range(len(parts)):
        cliqmatr = kronFor(symptop(rpaulis[parts[i][0]]))
        for j in range(1, len(parts[i])):
            cliqmatr = cliqmatr + kronFor(symptop(rpaulis[parts[i][j]]))
        out.append(cliqmatr)
    return out

def clique_trotterization(rpaulis, tset, toshow):
    ops = clique_matrices(rpaulis)
    oneset = []
    twoset = []
    threeset = []
    for time in tset:
        ideal = idealexp(ops, time)
        one = math.log(trdist(ideal, recurse(2, time, ops)), 10)
        two = math.log(trdist(ideal, recurse(4, time, ops)), 10)
        three = math.log(trdist(ideal, recurse(6, time, ops)), 10)
        oneset.append(one)
        twoset.append(two)
        threeset.append(three)
    print("usman array ; ", type(tset))
    output = [tset.tolist(), oneset, twoset, threeset]
    return output

def clique_savings_plot(npaulismax, nregs, toshow):
    xpts = []
    cliques = []
    estimate = []
    ratio = []
    physical_counts_nocliques_m22 = []
    physical_counts_cliques_m22 = []
    physical_counts_nocliques_m24 = []
    physical_counts_cliques_m24 = []
    physical_counts_nocliques_m26 = []
    physical_counts_cliques_m26 = []
    for i in range(10, npaulismax, int((npaulismax - 10) / 20)):
        cliques.append(len(partition(symptoadj(sympmatrix(randompaulis(i, nregs))))))
        xpts.append(i)
        estimate.append(i / math.log(i))
        ratio.append(cliques[-1] / estimate[-1])
        physical_counts_nocliques_m22.append(math.log(2 * i - 1, 10))
        physical_counts_cliques_m22.append(math.log(2 * cliques[-1] - 1, 10))
        physical_counts_nocliques_m24.append(math.log(math.pow(2 * i - 1, 5), 10))
        physical_counts_cliques_m24.append(math.log(math.pow(2 * cliques[-1] - 1, 5), 10))
        physical_counts_nocliques_m26.append(math.log(math.pow(2 * i - 1, 10), 10))
        physical_counts_cliques_m26.append(math.log(math.pow(2 * cliques[-1] - 1, 10), 10))
    output = [physical_counts_nocliques_m22,
              physical_counts_cliques_m22, physical_counts_nocliques_m24, physical_counts_cliques_m24, physical_counts_nocliques_m26, physical_counts_cliques_m26]
    return output

def streamnocheck(npaulis, nregs):
    paulis = randompaulis(npaulis, nregs)
    adjgraph = symptoadj(sympmatrix(paulis))
    cliques = recover_labels(partition(adjgraph), npaulis)
    ncliques = len(cliques)
    table = []
    for i in range(5):
        nocliques = math.pow(2 * npaulis - 1, 5 * (i + 1))
        withcliques = math.pow(2 * ncliques - 1, 5 * (i + 1))
        ratio = withcliques / nocliques
        nocliques = '{:0.3e}'.format(nocliques)
        withcliques = '{:0.3e}'.format(withcliques)
        ratio = '{:0.3e}'.format(ratio)
        table.append([nocliques, withcliques, ratio])
    output = clique_savings_plot(npaulis, nregs, False)
    result = {
        'length': len(cliques),
        'cliques_data': cliques,
        'table': table,
        'normal': output,
    }
    print("stream no check image")
    return result

def streamyescheck(npaulis, nregs):
    paulis = randompaulis(npaulis, nregs)
    adjgraph = symptoadj(sympmatrix(paulis))
    cliques = recover_labels(partition(adjgraph), npaulis)
    ncliques = len(cliques)
    table = []
    for i in range(5):
        nocliques = math.pow(2 * npaulis - 1, 5 * (i + 1))
        withcliques = math.pow(2 * ncliques - 1, 5 * (i + 1))
        ratio = withcliques / nocliques
        nocliques = '{:0.3e}'.format(nocliques)
        withcliques = '{:0.3e}'.format(withcliques)
        ratio = '{:0.3e}'.format(ratio)
        table.append([nocliques, withcliques, ratio])
    output = clique_savings_plot(npaulis, nregs, False)
    tset = numpy.arange(0.01, 0.1, 0.005)
    error = clique_trotterization(paulis, tset, False)
    result = {
        'length': len(cliques),
        'cliques_data': cliques,
        'table': table,
        'normal': output,
        'error': error,
    }
    print("stream yes check image")
    return result

