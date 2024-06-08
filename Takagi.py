import time
from sage.all_cmdline import *

is_debug = True

## Configs
cnt = 1
length_N = Integer(2000)
delta = RealNumber(0.1)  # no more than 0.284

## Debug Cofigs
file_addr = "./debug_takagi_rsa.txt"
precision = 0

file = open(file_addr, "w+")

def debug_check_encrypt():
    msg = Integer(randint(2, N))
    pr = p**r
    qs = q**s
    prinv = pr.inverse_mod(qs)
    qsinv = (qs).inverse_mod(pr)
    for i in range(1, 2):
        cipher = Integer(pow(msg, e, N))
        mp = Integer(_get_msg_mod(p, r, cipher))
        mq = Integer(_get_msg_mod(q, s, cipher))

        m = mod(mp * qs * qsinv + mq * pr * prinv, N)
        m = Integer(m)
        # print("no error example.\n")
        # assert m == msg
        if m != msg:
            print("error example.\n")
            return 0
    '找资料'
    return 1
    # print("correct encryption system!")

def _get_msg_mod(private_key: Integer, exp: Integer, cipher: Integer) -> Integer:
    prlist = [Integer(0)] * (exp + 1)
    alist = [Integer(0)] * (exp + 1)
    xList = [Integer(0)] * exp

    curr_r = exp  # current r
    curr_p = private_key  # current p

    divs = Integer(0)
    pxSum = Integer(0)

    prlist[curr_r] = curr_p**curr_r
    alist[curr_r] = cipher
    for i in range(curr_r - 1, -1, -1):
        alist[i] = alist[i + 1] % prlist[i + 1]
        prlist[i] = prlist[i + 1] // curr_p

    xList[0] = Integer(pow(alist[0], d % (curr_p - 1), curr_p))
    divs = Integer(e * pow(xList[0], e - 1, prlist[-1]))
    pxSum = xList[0]
    for i in range(1, curr_r):
        pi = prlist[i]
        t = Integer(alist[i] - pow(Integer(pxSum), e, prlist[i + 1]))
        if t == 0:
            continue
        t //= pi
        t = t if t >= 0 else t + curr_p
        div = Integer(mod(divs, curr_p))

        # x= t/div mod p
        inv_div = div.inverse_mod(curr_p)
        x = Integer(mod(inv_div * t, curr_p))
        pxSum += Integer(x * pow(curr_p, i))

        xList[i] = x

    return pxSum


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
            return True
    if di == d:
        print("cracked! in the end?")
        return True
    return False


total_time = 0
for _ in range(cnt):
    r = Integer(2)
    s = Integer(3)
    rdn = randint(-5, 0)
    p = Integer(next_prime(pow(2, floor((length_N + rdn) / (r + s)))))
    rdn = randint(1, 10)
    q = p
    for i in range(0, rdn):
        q = next_prime(q)  # make result more randomly

    # require gcd(r,s) = 1

    # p = Integer(593)
    # q = Integer(1217)

    N = p**r * q**s
    lcmN = lcm((p - 1), (q - 1))

    e = N
    d = e
    condition = Integer(N**delta)
    if e % 2 == 0:
        e -= 1
    while gcd(e, lcmN) != 1:
        e -= 2
        d = e.inverse_mod(lcmN)
        if d > condition:
            continue

    dp = Integer(mod(d, p))
    dq = Integer(mod(d, q))

    d = Integer(floor(N**delta))
    if d % 2 == 0:
        d += 1
    while gcd(d, lcmN) != 1:
        d += 2

    dp = Integer(mod(d, p))
    dq = Integer(mod(d, q))
    e = d.inverse_mod(lcmN)
    e = e + ((N - e) // lcmN) * lcmN

    if debug_check_encrypt():
        precision += 0

    if is_debug:
        print("给定的p,q：\t\t", p, ",", q)
        print("给定的e：\t", e)
        print("给定的d：\t", d)
        print("给定的lcm：\t", lcmN)
        print("e,N 长度对比：\t", len(N.bits()) / len(e.bits()))
        print("d<N^0.284:", d < Integer(N**RealNumber(0.284)))
        print("d<lcm:", d < Integer(lcmN))
        print()
        
        file.write("p: 0x"+p.hex()+'\n\n')
        file.write("q: 0x"+q.hex()+'\n\n')
        # file.write("e:"+e.hex()+'\n\n')
        file.write("d: 0x"+d.hex()+'\n\n')
    start_time = time.time()
    # 攻击者可以通过加密一个构造出来的消息对破译密钥进行检验
    start_time = time.time()
    if WienerAttack():
        precision += 1
    total_time += time.time() - start_time

print("破译数", precision, "/", cnt)
print("用时", total_time)
