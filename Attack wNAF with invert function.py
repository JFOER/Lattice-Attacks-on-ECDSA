from gmpy2 import invert
from fpylll import *
import time
import random
from fractions import Fraction

# Converting k into its window Non-Adjacent Form (wNAF) representation
# and deriving a perfect Double-Add-Invert chain.

def trans_wNAF(k, w):
    i = 0
    tmp = 0
    wNAF = []
    nz = []
    symbol = []
    while (k > 0):
        if k % 2 == 1:
            tmp = k % (2 ** (w + 1))
            if tmp >= 2 ** w:
                tmp = tmp - (2 ** (w + 1))
            wNAF.append(tmp)
            if tmp > 0:
                symbol.append("+")
            else:
                symbol.append("-")

            k = k - tmp
            nz.append(i)
        else:
            tmp = 0
            wNAF.append(tmp)
        k = k // 2
        i = i + 1
    return wNAF, nz, symbol


# Converting the Double-Add-Invert chain into leaked binary bits, as detailed in Section VI-B of the article.
def get_known_bits(k_nz, w):
    known = []
    if (k_nz[0] > 0):
        for i in range(0, k_nz[0]):
            known.append(i)
    for j in range(1, len(k_nz)):
        dis = k_nz[j] - k_nz[j - 1]
        for jj in range(w, dis):
            known.append(k_nz[j - 1] + jj)
        if (k_nz[j] <= 255):
            known.append(k_nz[j])
    return known


# Get the positions of unknown bits along with the bitsize and the sum of the known bits
# k: the nonce
# know: the leaked bits of nonce
# lam: the positions of unknown chunks
# mu: the bitsize of unknown chunks
# sum: the sum of known bits
# maxu: the maximum bitsize of unknown chunks(It is fixed referring to LSBs or MSBs)
def get_uk_sum(k, know):
    bin_k = bin(k)[2:]
    bin_rek = bin_k[::-1]
    unknow = []
    lam = []
    mu = []
    sum = 0
    maxu = 0
    for i in range(0, 256):
        if i not in know:
            unknow.append(i)
    i = 0
    while (i < len(unknow)):
        lam.append(unknow[i])
        count = 1
        while (i + count < len(unknow) and unknow[i] + count in unknow):
            count = count + 1
        mu.append(count)
        if (maxu < count):
            maxu = count
        i = i + count
    for j in know:
        if j < len(bin_rek):
            if (bin_rek[j] == "1"):
                sum = sum + 2 ** j
    return lam, mu, sum, maxu


# Randomly pick “u” signatures out of 1000 signatures for “num” times
def get_randomsig(u, num):
    numbers = range(1000)
    random_pairs = []
    while len(random_pairs) < num:
        rand_set = set()
        count = 0
        while (count < u):
            tmp = random.choice(numbers)
            while (tmp in rand_set):
                tmp = random.choice(numbers)
            rand_set.add(tmp)
            count = count + 1
        if rand_set not in random_pairs:
            random_pairs.append(rand_set)
    return random_pairs

# Merging method mentioned in Section IV-A
# lam: the positions of unknown chunks
# mu: the bitsize of unknown chunks
# limit: Set the limitation for loss data
# lam1: new lam after merging
# mu1: new mu after merging
# maxu: new maxu after merging
def merge(lam, mu, limit):
    lam1 = []
    mu1 = []
    maxu = 1
    loss1 = 0
    i = 0
    while (loss1 < limit and i <= len(lam) - 2):
        lam1.append(lam[i])
        flag = 0
        if (lam[i + 1] - mu[i] - lam[i] == 2):
            loss1 = loss1 + 1
            flag = 1
            tmp = i
            i = i + 1
            while (loss1 < limit and i <= len(lam) - 2 and lam[i + 1] - mu[i] - lam[i] == 1):
                loss1 = loss1 + 1
                i = i + 1
            mu1.append(lam[i] - lam[tmp] + mu[i])
            if (lam[i] - lam[tmp] + mu[i] > maxu):
                maxu = lam[i] - lam[tmp] + mu[i]
        if (flag == 0):
            mu1.append(mu[i])
            if (mu[i] > maxu):
                maxu = mu[i]
        i = i + 1
    while (i < len(mu)):
        lam1.append(lam[i])
        mu1.append(mu[i])
        i = i + 1
    return lam1, mu1, maxu


# Construct the lattice matrix
# gamma,T,tau,sigma,m,delta are the variable mentioned in paper(same notations)
# L is the dimension of the lattice which is T+1
# fac is the Z which defined in paper
# guesum is "correct guess"
def build_basis(gamma, u, T, L, m, tau, delta, sigma, k_mu, guesum, fac):
    A = IntegerMatrix(L, L)
    A[L - 1, L - 1] = 2 ** (m - 1)
    index = u - 2 + len(tau[0]) + 1
    pos = 0
    guess = []
    for i in range(0, u - 1):
        tmpx = (guesum[i]) % fac[i]
        guess.append(tmpx)
    print(f"The correct guess is {guess}")
    for i in range(0, u - 1):
        tmpmu = k_mu[i + 1][0]
        A[i, i] = fac[i] * q * delta * (2 ** (m - tmpmu))
        for j in range(0, len(tau[0])):
            A[u - 2 + j + 1, i] = tau[i][j] * delta * (2 ** (m - tmpmu)) + (2 ** (m - tmpmu)) * q
            A[u - 2 + j + 1, u - 2 + j + 1] = 2 ** (m - k_mu[0][j])
            A[L - 1, u - 2 + j + 1] = 2 ** (m-1)
        for tmp in range(0, len(sigma[i])):
            A[index, i] = delta * (2 ** (m - tmpmu)) * sigma[i][tmp]
            A[index, index] = 2 ** (m - k_mu[i + 1][pos + 1])
            A[L - 1, index] = 2 ** (m-1)
            for xindex in range(0, u - 1):
                tmpmmu = k_mu[xindex + 1][0]
                A[index, xindex] = A[index, xindex] + (2 ** (m - tmpmmu)) * q
            index = index + 1
            pos = pos + 1
        pos = 0
        A[L - 1, i] = delta * (2 ** (m - tmpmu)) * gamma[i] + 2 ** (m-1)  + q * (2 ** (m - tmpmu)) * guess[i]
    return A


#index: There are 100 groups of data.7z, pick one (e.g., 29)
index=29
#q: Modulus of elliptic curve(Curve secp256k1)
q = 115792089237316195423570985008687907852837564279074904382605163141518161494337


filename_k = "/home/jfoer/data/data" + str(index) + "/k.txt"
filename_s = "/home/jfoer/data/data" + str(index) + "/s.txt"
filename_r = "/home/jfoer/data/data" + str(index) + "/r.txt"
filename_hash = "/home/jfoer/data/data" + str(index) + "/hash.txt"
filename_pri = "/home/jfoer/data/data" + str(index) + "/pri_key.txt"

# k,hash,r,s: signature data, 1000 groups in total
k = []
hash = []
r = []
s = []
#pri_key : private key
pri_key = 0

k_know = []
k_bar = []
k_lam = []
k_mu = []
k_m = []
k_wNAF=[]
k_nz=[]
k_symbol=[]
w=3

with open(filename_pri, 'r') as file:
    for line in file:
        line = line.strip()
        pri_key = int(line, 16)

with open(filename_k, 'r') as file:
    i = -1
    for line in file:
        i = i + 1
        line = line.strip()
        k_decimal = int(line, 16)
        k.append(k_decimal)

with open(filename_r, 'r') as file:
    j = -1
    for line in file:
        j = j + 1
        line = line.strip()
        r_decimal = int(line, 16)
        r.append(r_decimal)

with open(filename_s, 'r') as file:
    j = -1
    for line in file:
        j = j + 1
        line = line.strip()
        s_decimal = int(line, 16)
        s.append(s_decimal)

with open(filename_hash, 'r') as file:
    j = -1
    for line in file:
        j = j + 1
        line = line.strip()
        hash_decimal = int(line, 16)
        hash.append(hash_decimal)

for i in range(0,1000):
        tmp1,tmp2,tmp3=trans_wNAF(k[i],w)
        k_wNAF.append(tmp1)
        k_nz.append(tmp2)
        k_symbol.append(tmp3)
for i in range(0,1000):
    tmp1=get_known_bits(k_nz[i],w)
    k_know.append(tmp1)
    tmp2,tmp3,tmp4,tmp5=get_uk_sum(k[i],tmp1)
    k_lam.append(tmp2)
    k_mu.append(tmp3)
    k_bar.append(tmp4)
    k_m.append(tmp5)

# u: the number of signatures
u = 3
# num: the number of the experiments
num = 200
random_pairs = get_randomsig(u, num)
# suc: The number of successful times
suc = 0
sum_time = 0
for cnt in range(0, num):
    sig_pos = random_pairs[cnt]
    t_k = []
    t_k_bar = []
    t_k_lam = []
    t_k_mu = []
    t_k_m = []
    t_hash = []
    t_r = []
    t_s = []
    t_know = []
    for item in sig_pos:
        print(f"Use the signature with index {item}")
        t_k.append(k[item])
        t_know.append(k_know[item])
        t_k_bar.append(k_bar[item])
        t_k_lam.append(k_lam[item])
        t_k_mu.append(k_mu[item])
        t_k_m.append(k_m[item])
        t_hash.append(hash[item])
        t_r.append(r[item])
        t_s.append(s[item])
    mr_lam = []
    mr_mu = []
    #limit: Set the limitation for loss data
    limit = 20

    m = 1
    for i in range(0, u):
        if (t_k_m[i] > m):
            m = t_k_m[i]
    for i in range(0, u):
        tmp1, tmp2, tmp3 = merge(t_k_lam[i], t_k_mu[i], limit)
        if tmp3 > m:
            m = tmp3
        mr_lam.append(tmp1)
        mr_mu.append(tmp2)
    gamma = []
    tau = []
    sigma = []
    # caculate the gamma,tau,sigma
    for i in range(1, u):
        tmprev = t_s[i] * t_r[0] * 2 ** (mr_lam[i][0])
        rev = int(invert(tmprev,q))
        tmp3 = []
        for kk in mr_lam[i]:
            tmp3.append((-1) * rev * t_s[i] * t_r[0] * 2 ** kk)
        sigma.append(tmp3[1:])
        tmp1 = rev * ((t_r[0] * (t_s[i] * t_k_bar[i] - t_hash[i])) - (t_r[i] * (t_s[0] * t_k_bar[0] - t_hash[0])))
        gamma.append(tmp1)
        tmp2 = []
        for j in mr_lam[0]:
            tmp2.append(rev * t_s[0] * t_r[i] * (2 ** j))
        tau.append(tmp2)

    T = 0
    for i in range(0, u):
        T = T + len(mr_lam[i])
    # L: The dimension of the lattice
    L = T + 1
    print(f"The dimension of the lattice is {L}")
    #For "correct guess"
    dv = []
    for i in range(0, u):
        tmptm = bin(t_k[i])[2:]
        tmp = tmptm[::-1]
        tmpx = []
        for j in range(0, len(mr_mu[i])):
            tmp1 = 0
            for zx in range(0, mr_mu[i][j]):
                if (mr_lam[i][j] + zx < len(tmp) and tmp[mr_lam[i][j] + zx] == "1" and mr_lam[i][j] + zx not in
                        t_know[i]):
                    tmp1 = tmp1 + 2 ** (zx)
            tmpx.append(tmp1)
        dv.append(tmpx)
    dsum = []
    for i in range(0, u):
        tmpdsum = 0
        if (i == 0):
            tmpdsum = tmpdsum + dv[0][0]
        for j in range(1, len(dv[i])):
            tmpdsum = tmpdsum + dv[i][j]

        dsum.append(tmpdsum)
    modt = []
    for i in range(1, u):
        tmprevd = t_s[i] * t_r[0] * 2 ** (mr_lam[i][0])
        revd = int(invert(tmprevd,q))
        xxx = ((tmprevd * revd * dv[i][0]) // q)
        tmpt = (revd * (t_r[i] * t_s[0] * t_k[0] - t_r[0] * t_s[i] * t_k[i] + t_r[0] * t_hash[i] - t_r[i] * t_hash[
            0] + tmprevd * dv[i][0])) // q
        modt.append(tmpt)
    guesum = []
    dexp0 = 0
    for i in range(1, u):
        dexp0 = dexp0 + dsum[i]
    for i in range(0, u - 1):
        tmpgs = dsum[0] + dexp0 + modt[i]
        guesum.append(tmpgs)

    delta = 1
    #fac: It's Z mentioned in paper
    #For different leakage models, it is recommended to use the Z mentioned in the paper.
    #Or choose other values, and the higher the value, the higher the success rate.
    fac = [2 ** 9]
    tmplen = len(fac)
    for i in range(tmplen, u - 1):
        fac.append(1)
    lattice = build_basis(gamma, u, T, L, m, tau, delta, sigma, mr_mu, guesum, fac)
    begin_time = time.time()
    param = BKZ.Param(block_size = 45, strategies = BKZ.DEFAULT_STRATEGY )
    l_bkz = BKZ.reduction(lattice, param)
    end_time = time.time() - begin_time
    sum_time = sum_time + end_time
    flag = 0
    r0inv = int(invert(t_r[0],q))
    beta0 = (t_s[0] * t_k_bar[0] - t_hash[0]) % q
    for i in range(0, L):
        if (l_bkz[i][L - 1] == 2 ** (m - 1)):
            tmp_sum1 = 0
            for index2 in range(0, len(mr_lam[0])):
                tmpfrac = Fraction(2 ** mr_mu[0][index2], 2 ** m)
                tmp_sum1 = tmp_sum1 + (2 ** (mr_lam[0][index2]) * t_s[0]) * (
                            (-l_bkz[i][u - 1 + index2] + 2 ** (m - 1)) * tmpfrac)

            alpha = ((tmp_sum1 + beta0) * r0inv) % q

            if (alpha == pri_key):
                flag = 1
                break
        if (l_bkz[i][L - 1] == -2 ** (m - 1)):
            tmp_sum1 = 0
            for index2 in range(0, len(mr_lam[0])):
                tmpfrac = Fraction(2 ** mr_mu[0][index2], 2 ** m)
                tmp_sum1 = tmp_sum1 + (2 ** (mr_lam[0][index2]) * t_s[0]) * (
                            (l_bkz[i][u - 1 + index2] + 2 ** (m - 1)) * tmpfrac)
            alpha = ((tmp_sum1 + beta0) * r0inv) % q
            if (alpha == pri_key):
                flag = 1
                break
    if (flag == 1):
        print(f"delta is {delta}")
        print(f"The Duration of this experiment {end_time}s")
        print(f"The Average Duration till this time {sum_time / (cnt + 1)}s")
        print("correct!")
        suc = suc + 1
        print(f"This is the data from index{index}, Test number: {cnt + 1}, Current number of successful attempts: {suc},total: {num}")
    if (flag == 0):
        print(f"delta is {delta}")
        print(f"The Duration of this experiment {end_time}s")
        print(f"The Average Duration till this time {sum_time / (cnt + 1)}s")
        print("wrong answer!")
        print(f"This is the data from index{index}, Test number: {cnt + 1}, Current number of successful attempts: {suc},total: {num}")
print(f"The success rate is {float(suc / (num//100))}%")
#Please note that the formula for calculating the average duration mentioned in the paper is: Z / (32 (threads) * 2) * Average duration.
