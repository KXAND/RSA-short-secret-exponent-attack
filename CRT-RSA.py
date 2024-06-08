import random
import time
from sage.all_cmdline import *
from math import log2

is_debug = True

## Configs
cnt = 1
length_N = Integer(4096)
delta = RealNumber(0.24)  # no more than 0.25

## Debug Cofigs
file_addr = "debug_crt_rsa.txt"
precision = 0

file = open(file_addr, "w+")


def debug_check_encrypt(N, e, dp, dq):
    msg = Integer(randint(2, N))
    pinv = p.inverse_mod(q)
    qinv = q.inverse_mod(p)
    while i in range(1, 50):
        cipher = pow(msg, e, N)
        cp = Integer(mod(cipher, p))
        cq = Integer(mod(cipher, q))
        mp = Integer(pow(cp, dp, p))
        mq = Integer(pow(cq, dq, q))

        m = mod(mp * q * qinv + mq * p * pinv, N)
        m = Integer(m)
        print("no error example.\n")
        assert m == msg
    print("correct encryption system!")


"""
Returns:
return denominators and k-th continued fraction sum
"""


def rational_to_continued_fraction(e, N, coefficients: list = None):
    convergents = [(1, 0)]  # helper for convergent calculation

    q = e // N
    if coefficients:
        coefficients.append(q)  # 系数
    convergents.append((q, 1))  # 收敛项

    convpp = (1, 0)  # previous's previous one of convergents
    while q * N != e:
        temp = N
        N = e - q * N
        e = temp
        q = e // N  # q_i = q_{i-1}的余数取倒再整除
        if coefficients:
            coefficients.append(q)
        cn = q * convergents[-1][0] + convpp[0]
        cd = q * convergents[-1][1] + convpp[1]
        convpp = convergents[-1]
        convergents.append((cn, cd))
    return convergents


def WienerAttack():
    convs = rational_to_continued_fraction(e, N, None)
    i = 0
    for ki, di in convs:
        i += 1
        if not ki:
            continue

        phiN = (e * di - 1) // ki
        if phiN == 0:
            continue

        # pqSum = N - phiN + 1  # p+q
        # if pqSum & 1:
        #     continue

        # pqDiff = pqSum * pqSum - 4 * N  # p-q
        # if pqDiff < 0:
        #     continue
        # pqDiff = floor(sqrt(pqDiff))

        # if pqDiff**2 != pqSum**2 - 4 * N:  # perfect square
        #     continue

        if di == d:
            print("cracked!")
            print("find d at ", i - 2, " of ", len(convs) - 2)
            file.write("find d at " + str(i - 2) + " of " + str(len(convs) - 2)+'\n')
            return True
    if di == d:
        print("cracked! in the end?")
        return True
    return False


total_time = 0
for _ in range(cnt):
    rdn = randint(-5, 0)
    p = Integer(next_prime(pow(2, floor((length_N + rdn) / 2))))
    rdn = randint(1, 10)
    q = p
    for i in range(0, rdn):
        q = next_prime(q)  # make result more randomly

    N = p * q
    phiN = (p - 1) * (q - 1)
    d = Integer(floor(N**delta))
    if d % 2 == 0:
        d += 1
    while gcd(d, phiN) != 1:
        d += 2
    dp = Integer(mod(d, p))
    dq = Integer(mod(d, q))
    e = d.inverse_mod(phiN)

    if is_debug:
        print("给定的p,q：\t\t", p, ",", q)
        print("d>p?", d > p)
        file.write("p: 0x"+p.hex()+'\n')
        file.write("q: 0x"+q.hex()+'\n')
        # file.write("e: 0x"+e.hex()+'\n')
        file.write("d: 0x"+d.hex()+'\n')
    start_time = time.time()
    # 攻击者可以通过加密一个构造出来的消息对破译密钥进行检验
    start_time = time.time()
    if WienerAttack():
        precision += 1
    total_time += time.time() - start_time

print("破译数", precision, "/", cnt)
print("用时", total_time)
