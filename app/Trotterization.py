import numpy
import math
import scipy
import matplotlib
import scipy.linalg
import matplotlib.pyplot as pyplot
import cmath
from datetime import datetime

##finally installed correctly


# returns the trace distance between operators A and B
def trdist(A, B):
    U, D, V = numpy.linalg.svd(numpy.subtract(A, B))
    return numpy.sum(D)


def thirdlevel(A, B, t):
    s = 1. / (2 - math.pow(2, 1. / 3))
    # print(s)
    rollingproduct = scipy.linalg.expm((s / 2.) * t * A)
    rollingproduct = numpy.matmul(scipy.linalg.expm(s * t * B), rollingproduct)
    rollingproduct = numpy.matmul(scipy.linalg.expm((1 - s) / 2. * t * A), rollingproduct)
    rollingproduct = numpy.matmul(scipy.linalg.expm((1 - 2 * s) * t * B), rollingproduct)
    rollingproduct = numpy.matmul(scipy.linalg.expm((1 - s) / 2. * t * A), rollingproduct)
    rollingproduct = numpy.matmul(scipy.linalg.expm(s * t * B), rollingproduct)
    rollingproduct = numpy.matmul(scipy.linalg.expm(s / 2. * t * A), rollingproduct)
    return rollingproduct


def firsttest(A, B, t):
    ideal = scipy.linalg.expm(t * numpy.add(A, B))
    badapprox = numpy.matmul(scipy.linalg.expm(t * A), scipy.linalg.expm(t * B))
    betterapprox = numpy.matmul(numpy.matmul(scipy.linalg.expm((t / 2.) * A), scipy.linalg.expm(t * B)),
                                scipy.linalg.expm((t / 2.) * A))
    one = trdist(ideal, badapprox)
    two = trdist(ideal, betterapprox)
    three = trdist(ideal, thirdlevel(A, B, t))
    # print(one)
    # print(two)
    ##add third order live, from suzuki paper
    ptsone = badapprox
    ptstwo = betterapprox
    ptsthree = thirdlevel(A, B, t)
    return [one, two, three, ptsone, ptstwo, ptsthree]


def firstplot(A, B, tset):
    oneset = []
    twoset = []
    threeset = []
    ptsone = []
    ptstwo = []
    ptsthree = []
    example = numpy.array([1, 1])
    for time in tset:
        temp = firsttest(A, B, time)
        oneset.append(temp[0])
        twoset.append(temp[1])
        threeset.append(temp[2])
        ptsone.append(temp[3].dot(example))
        ptstwo.append(temp[4].dot(example))
        ptsthree.append(temp[5].dot(example))
    # print(oneset)
    # print(twoset)
    pyplot.plot(tset, oneset, color='tab:blue')
    pyplot.plot(tset, twoset, color='tab:orange')
    pyplot.plot(tset, threeset, color='tab:red')
    #    pyplot.plot(tset,ptsone,"b*",tset,ptstwo,"rp",tset,ptsthree,"c^")
    pyplot.show()
    return


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
    ##make it match the standard convention
    paulis = str(paulis[::-1])
    res = numpy.kron(matChoose(paulis[1]), matChoose(paulis[0]))

    for i in range(2, len(paulis)):
        ##flipped order for same reason
        res = numpy.kron(matChoose(paulis[i]), res)

    return res


# reverse product
# scipy.linalg.expm stuff

# converts paulis to matrices
def paulistomatrices(paulis):
    matr = []
    for i in range(len(paulis)):
        matr.append(kronFor(paulis[i]))
    return matr


# has reverse prod
# change paulis to ops and feed in results from kronFor
##def TrotMoreArgs(paulis, x):
##    res = numpy.identity(pow(2, len(paulis[0])))
##    for i in range(len(paulis)):
##        mat = kronFor(paulis[i])
##        mat = mat * x * 0.5
##        fin = scipy.linalg.expm(mat)
##        res = res.dot(fin)
##
##    for i in range(len(paulis)-1, -1, -1):
##        mat = kronFor(paulis[i])
##        mat = mat * x * 0.5
##        fin = scipy.linalg.expm(mat)
##        res = res.dot(fin)
##
##    return res

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


# needs updating to handle more/general args
##def idealexp(paulis,x):
##    res= numpy.matrix(kronFor(paulis[0]))
##    for i in range(1,len(paulis)):
##        res=numpy.add(res,numpy.matrix(kronFor(paulis[i])))
##    return scipy.linalg.expm(x*res)

def idealexp(ops, x):
    res = numpy.matrix(ops[0])
    for i in range(1, len(ops)):
        res = numpy.add(res, numpy.matrix(ops[i]))
    return scipy.linalg.expm(x * res)


# print(kronFor("ixz"))

def kronAll(paulis, index):
    if index == len(paulis):
        return
    if index == 1:
        res = numpy.kron(matChoose(paulis[index - 1]), matChoose(paulis[index]))
        kronAll(paulis, index + 1)

    res = numpy.kron(res, matChoose(paulis[index]))
    kronAll(paulis, index + 1)


###these functions need altering so that I can put in arbitrary matrices, not just paulis
###I want to put in sums of Paulis
##replace paulis with ops?
##def recurse(m2, x, paulis):
##    if(m2 <= 2):
##        return TrotMoreArgs(paulis, x)
##    p = pow((4 - pow(4, (1. / (m2 - 1)))), -1)
##    fin = (numpy.linalg.matrix_power(recurse((m2-2), (p*x), paulis),2).dot(recurse((m2-2),((1-(4*p))*x), paulis))).dot(numpy.linalg.matrix_power(recurse((m2-2),(p*x),paulis),2))
##    return fin

def recurse(m2, x, ops):
    if (m2 <= 2):
        return TrotMoreArgs(ops, x)
    p = pow((4 - pow(4, (1. / (m2 - 1)))), -1)
    fin = (numpy.linalg.matrix_power(recurse((m2 - 2), (p * x), ops), 2).dot(
        recurse((m2 - 2), ((1 - (4 * p)) * x), ops))).dot(numpy.linalg.matrix_power(recurse((m2 - 2), (p * x), ops), 2))
    return fin


##def newplot(paulis,tset):
##    oneset = []
##    twoset = []
##    threeset = []
##    for time in tset:
##        ideal=idealexp(paulis,time)
##        one=math.log(trdist(ideal,recurse(2,time,paulis)))
##        two=math.log(trdist(ideal,recurse(4,time,paulis)))
##        three=math.log(trdist(ideal,recurse(6,time,paulis)))
##        oneset.append(one)
##        twoset.append(two)
##        threeset.append(three)
##    pyplot.plot(tset, oneset,color ='tab:blue')
##    pyplot.plot(tset, twoset,color ='tab:orange')
##    pyplot.plot(tset, threeset,color ='tab:red')
##    pyplot.show()
##    return

def newplot(ops, tset):
    oneset = []
    twoset = []
    threeset = []
    for time in tset:
        ideal = idealexp(ops, time)
        one = math.log(trdist(ideal, recurse(2, time, ops)))
        two = math.log(trdist(ideal, recurse(4, time, ops)))
        three = math.log(trdist(ideal, recurse(6, time, ops)))
        oneset.append(one)
        twoset.append(two)
        threeset.append(three)
    pyplot.plot(tset, oneset, color='tab:blue')
    pyplot.plot(tset, twoset, color='tab:orange')
    pyplot.plot(tset, threeset, color='tab:red')
    pyplot.show()
    return


####note: thirdlevel gives S_3 not S_3^*!


# computes comm of two paulis
def comm(p1, p2):
    comm = 0
    n = int(len(p1) / 2)
    for i in range(n):
        comm = comm + (p1[i] * p2[i + n] - p1[i + n] * p2[i])
    return comm % 2


# returns the full symp prod matrix
def sympmatrix(paulis):
    m = numpy.zeros([len(paulis), len(paulis)])
    for i in range(len(m)):
        for j in range(len(m)):
            m[i, j] = comm(paulis[i], paulis[j])
    return m


# returns the adjacency matrix for the graph from the paulis
def symptoadj(symp):
    for i in range(len(symp)):
        for j in range(len(symp)):
            if (not (i == j)):
                symp[i, j] = 1 - symp[i, j]
    return symp.tolist()


import random


# takes in an adjacency matrix
# returns one maximal clique from the graph
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


# compresses adj to columns and rows not including set vrem
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
        # vertex ids are messed up from remove, so use recover_labels
    return parts


##returns v, non-identical paulis on n registers
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


####have to return clique counts
###script to sum clique members
###take clique indices, transform to paulis then to matrix then sum

# takes random paulis (rpaulis) and returns the matrices associated to the cliques
def clique_matrices(rpaulis):
    out = []
    adjmat = symptoadj(sympmatrix(rpaulis))
    parts = recover_labels(partition(adjmat), len(rpaulis))
    # print(len(parts))
    for i in range(len(parts)):
        cliqmatr = kronFor(symptop(rpaulis[parts[i][0]]))
        for j in range(1, len(parts[i])):
            cliqmatr = cliqmatr + kronFor(symptop(rpaulis[parts[i][j]]))
        out.append(cliqmatr)
    return out


###add new function that feeds in random paulis and generates the plots
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
    # pyplot.plot(tset, oneset, color='tab:blue', linewidth=4)
    # pyplot.plot(tset, twoset, color='tab:orange', linewidth=4)
    # pyplot.plot(tset, threeset, color='tab:red', linewidth=4)
    # pyplot.rcParams.update({'font.size':18})
    # pyplot.rc('xtick',labelsize=48)
    # pyplot.xlabel('time', fontsize=24)
    # pyplot.ylabel('log_10 trace distance', fontsize=24)
    # pyplot.yticks(fontsize=24)
    # pyplot.xticks(fontsize=24)
    # if (toshow):
    #     pyplot.show()
    print("usman array ; ", type(tset))
    output = [tset.tolist(), oneset, twoset, threeset]
    return output


# m2 is 2*m
def physical_ops_count(rpaulis, m2):
    noclique = 2 * len(rpaulis) - 1
    clique = 2 * len(partition(symptoadj(sympmatrix(rpaulis)))) - 1
    for i in range(int(m2 / 2)):
        noclique = math.pow(noclique, 5)
        clique = math.pow(clique, 5)
    return [noclique, clique]


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
    #    pyplot.plot(xpts,cliques,linewidth=4)
    #    pyplot.plot(xpts,estimate,linewidth=4)
    # pyplot.plot(xpts,ratio)
    # pyplot.plot(xpts, physical_counts_nocliques_m22, "c-*", linewidth=2)
    # pyplot.plot(xpts, physical_counts_cliques_m22, "k-*", linewidth=2)
    # pyplot.plot(xpts, physical_counts_nocliques_m24, "c-p", linewidth=2)
    # pyplot.plot(xpts, physical_counts_cliques_m24, "k-p", linewidth=2)
    # pyplot.plot(xpts, physical_counts_nocliques_m26, "c-^", linewidth=2)
    # pyplot.plot(xpts, physical_counts_cliques_m26, "k-^", linewidth=2)
    # pyplot.xlabel('Number of Paulis (6 qubits)', fontsize=24)
    # pyplot.ylabel('log_10 of Number of Physical Operations', fontsize=24)
    # pyplot.yticks(fontsize=24)
    # pyplot.xticks(fontsize=24)
    # if (toshow):
    #     pyplot.show()
    output = [xpts, physical_counts_nocliques_m22,
              physical_counts_cliques_m22, physical_counts_nocliques_m24, physical_counts_cliques_m24, physical_counts_nocliques_m26, physical_counts_cliques_m26]
    return output


def streamnocheck(npaulis, nregs):
    # if (npaulis > 1000):
    #     return "too large a set was requested"
    # if (npaulis >= math.pow(2, 2 * nregs)):
    #     return "too many operators for the number of registers"
    paulis = randompaulis(npaulis, nregs)
    adjgraph = symptoadj(sympmatrix(paulis))
    cliques = recover_labels(partition(adjgraph), npaulis)
    # print(len(cliques))
    # print(cliques)
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
    # figure automatically says 6 regs, can we have it use the given value?
    # savings_plot_nocheck = f"/media/graph/savings_plot_nocheck-{datetime.now().date()}-{datetime.now().timestamp()}.png"
    # pyplot.savefig('.'+savings_plot_nocheck)
    # pyplot.clf()
    # pyplot.close('all')
    result = {
        'length': len(cliques),
        'cliques_data': cliques,
        'table': table,
        'normal': output,
    }
    print("stream no check image")
    return result


def streamyescheck(npaulis, nregs):
    # if (npaulis > 100):
    #     return "too large a set was requested"
    # if (npaulis >= math.pow(2, 2 * nregs)):
    #     return "too many operators for the number of registers"
    paulis = randompaulis(npaulis, nregs)
    adjgraph = symptoadj(sympmatrix(paulis))
    cliques = recover_labels(partition(adjgraph), npaulis)
    # print(len(cliques))
    # print(cliques)
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
    # print(table)
    normal = clique_savings_plot(npaulis, nregs, False)
    # figure automatically says 6 regs, can we have it use the given value?
    # savings_plot_yes_check = f"/media/graph/savings_plot_yes_check-{datetime.now().date()}-{datetime.now().timestamp()}.png"
    # savings_plot_error_check = f"/media/graph/savings_plot_error_check-{datetime.now().date()}-{datetime.now().timestamp()}.png"
    # pyplot.savefig('.'+savings_plot_yes_check)
    tset = numpy.arange(0.01, 0.1, 0.005)
    # pyplot.clf()
    # pyplot.close('all')
    error = clique_trotterization(paulis, tset, False)
    # pyplot.savefig('.'+savings_plot_error_check)
    # pyplot.clf()
    # pyplot.close('all')
    result = {
        'length': len(cliques),
        'cliques_data': cliques,
        'table': table,
        'normal': normal,
        'error': error,
    }
    print("stream yes check image")
    return result

# npaulis = 50
# registers = 3
# graphy = randompaulis(npaulis, registers)
# sympmat = sympmatrix(graphy)
# print(graphy)
# print(sympmat)
# print(symptoadj(sympmat))
##exampleadj=symptoadj(sympmat)


# exampleadj=[[0, 0, 1, 1, 0, 1],
# [0, 0, 1, 0, 0, 0],
# [1, 1, 0, 1, 1, 1],
# [1, 0, 1, 0, 1, 0],
# [0, 0, 1, 1, 0, 1],
# [1, 0, 1, 0, 1, 0]]


# exampleadj=[[0,1,1,0,1,0],
#            [1,0,1,1,0,1],
#            [1,1,0,1,0,1],
#            [0,1,1,0,1,1],
#            [1,0,0,1,0,0],
#            [0,1,1,1,0,0]]

# exampleadj=[[0]]

# temp=find_clique(exampleadj)
# print(temp)
# print(remove_verts(exampleadj,temp))

##cliquess=recover_labels(partition(exampleadj),npaulis)
##print("Number of cliques:")
##print(len(cliquess))
##print("Clique members")
##print(cliquess)
##print("Clique operators")
##print(clique_matrices(graphy))


# trotter method done recursively

# print(numpy.kron([[1,0],[0,1]],[[0,1],[1,0]]))
# print(trdist([[1,0],[0,-1]],[[0,0],[0,0]]))
# print(firsttest(numpy.array([[1.,0.],[0.,-1.]]),numpy.array([[0.,1.],[1.,0.]]),1.))

# print(cmath.exp(complex(0, 1) * cmath.pi))

# ts = numpy.arange(0.01, 0.1, 0.005)
###firstplot(complex(0, 1) * numpy.array([[1., 0.], [0., -1.]]), complex(0, 1) * numpy.array([[0., 1.], [1., 0.]]), ts)


###note: if these commute, the error is effectively 0--as we would expect
# arr = ["xx", "zx"]
# arr = ["z","x"]
# print(TrotMoreArgs(arr, 0.5))

####newplot(paulistomatrices(arr),ts)


# example_paulis=randompaulis(100,4)
# clique_trotterization(example_paulis,ts)
# print(physical_ops_count(example_paulis,6))

###clique_savings_plot(1000,6)

# streamyescheck(100, 4)
