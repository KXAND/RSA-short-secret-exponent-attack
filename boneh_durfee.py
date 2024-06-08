import time
from sage.all_cmdline import *

""" 
Returns:
* 0,0   if fail
* x0,y0 the solutions
"""


def boneh_durfee(pol, modulus, mm, tt, XX, YY, file):
    debug = True
    # substitution (Herrman and May)
    PR = PolynomialRing(
        ZZ,
        names=(
            "u",
            "x",
            "y",
        ),
    )
    (
        u,
        x,
        y,
    ) = PR._first_ngens(3)
    Q = PR.quotient(x * y + Integer(1) - u)  # u = xy + 1
    polZ = Q(pol).lift()

    UU = XX * YY + Integer(1)

    # x-shifts
    gg = []
    for kk in range(mm + Integer(1)):
        for ii in range(mm - kk + Integer(1)):
            xshift = x**ii * modulus ** (mm - kk) * polZ(u, x, y) ** kk
            gg.append(xshift)
    gg.sort()

    # x-shifts list of monomials
    monomials = []
    for polynomial in gg:
        for monomial in polynomial.monomials():
            if monomial not in monomials:
                monomials.append(monomial)
    monomials.sort()

    # y-shifts selected by Herrman and May
    for jj in range(Integer(1), tt + Integer(1)):
        for kk in range(floor(mm / tt) * jj, mm + Integer(1)):
            yshift = y**jj * polZ(u, x, y) ** kk * modulus ** (mm - kk)
            yshift = Q(yshift).lift()
            gg.append(yshift)

    # y-shifts list of monomials
    for jj in range(Integer(1), tt + Integer(1)):
        for kk in range(floor(mm / tt) * jj, mm + Integer(1)):
            monomials.append(u**kk * y**jj)
    # if debug:
        # file.write("=== parameters ===\n")
        # file.write("pol" + str(pol) + "\n")
        # file.write("e" + str(modulus) + "\n")
        # file.write("m" + str(mm) + "\n")
        # file.write("t" + str(tt) + "\n")
        # file.write("X" + str(XX) + "\n")
        # file.write("Y" + str(YY) + "\n")
        # file.write("=== polynomials ===\n")
        # for polynomial in gg:
        #     file.write(str(polynomial) + "\n")

    # construct lattice B
    nn = len(monomials)
    BB = Matrix(ZZ, nn)
    for ii in range(nn):
        BB[ii, Integer(0)] = gg[ii](Integer(0), Integer(0), Integer(0))
        for jj in range(Integer(1), ii + Integer(1)):
            if monomials[jj] in gg[ii].monomials():
                BB[ii, jj] = gg[ii].monomial_coefficient(monomials[jj]) * monomials[jj](
                    UU, XX, YY
                )

    BB = BB.LLL()

    if debug:
        print("LLL is done!")

    # transform vector i & j -> polynomials 1 & 2
    if debug:
        print("looking for independent vectors in lattice")
    found_polynomials = False

    for pol1_idx in range(nn - Integer(1)):
        for pol2_idx in range(pol1_idx + Integer(1), nn):
            # for i and j, create the two polynomials
            PR = PolynomialRing(
                ZZ,
                names=(
                    "w",
                    "z",
                ),
            )
            (
                w,
                z,
            ) = PR._first_ngens(2)
            pol1 = pol2 = Integer(0)
            for jj in range(nn):
                pol1 += (
                    monomials[jj](w * z + Integer(1), w, z)
                    * BB[pol1_idx, jj]
                    / monomials[jj](UU, XX, YY)
                )
                pol2 += (
                    monomials[jj](w * z + Integer(1), w, z)
                    * BB[pol2_idx, jj]
                    / monomials[jj](UU, XX, YY)
                )

            # resultant
            PR = PolynomialRing(ZZ, names=("q",))
            (q,) = PR._first_ngens(1)
            rr = pol1.resultant(pol2)

            if rr.is_zero() or rr.monomials() == [Integer(1)]:
                continue
            else:
                print("using vectors", pol1_idx, "and", pol2_idx)
                found_polynomials = True
                break
        if found_polynomials:
            break

    if not found_polynomials:
        print("no independant vectors could be found...")
        return 0, 0

    rr = rr(q, q)

    # solutions
    soly = rr.roots()

    if len(soly) == 0:
        print("delta is too small")
        return 0, 0

    # resultant
    soly = soly[0][0]
    ss = pol1(q, soly)
    solx = ss.roots()[0][0]

    return solx, soly


def crt_rsa_logic(file, length_N, length_d, delta, m):
    p = Integer(next_prime(2 ** int(round(length_N / 2))))
    # q = Integer(next_prime(round(pi.n() * p)))
    q = Integer(next_prime(next_prime(p)))
    N = Integer(p * q)
    phiN = (p - 1) * (q - 1)

    d = Integer(int(pow(N, length_d)))
    if d % 2 == 0:
        d += 1
    while gcd(d, phiN) != 1:
        d += 2
    e = d.inverse_mod(phiN)

    t = floor((1 - 2 * delta) * m)  # optimization from Herrmann and May
    Y = floor(pow(N, RealNumber(0.5)))
    X = floor(2 * floor(N**delta))
    A = int((N + 1) / 2)

    # s & k in k(A+s) = 1 mod(e)
    y0 = (-p - q) / 2  # s
    x0 = (1 - e * d) / (A + y0)  # k

    P = PolynomialRing(
        ZZ,
        names=(
            "x",
            "y",
        ),
    )
    (
        x,
        y,
    ) = P._first_ngens(2)
    pol = x * (A + y) + 1

    print("=== checking values ===")
    print("* |y| < Y:", abs(y0) < Y)
    print("* |x| < X:", abs(x0) < X)
    print("* t:", t)
    print("* d < N**0.292?:", d < N ** (0.292))
    print("* size of d:", int(floor(log(d)) / log(2)))

    # file
    if debug:
        file.write("p: 0x" + p.hex() + "\n\n")
        file.write("q: 0x" + q.hex() + "\n\n")
        file.write("d: 0x" + d.hex() + "\n\n")
        file.write("m: 0x" + m.hex() + "\n\n")
        file.write("t: 0x" + t.hex() + "\n\n")
        # file.write("=== checking values ===\n")
        # file.write("* p:" + str(p) + "\n")
        # file.write("* q:" + str(q) + "\n")
        # file.write("* phiN:" + str(phiN) + "\n")
        # file.write("* N:" + str(N) + "\n")
        # file.write("* e:" + str(e) + "\n")
        # file.write("* d:" + str(d) + "\n")
        # file.write("=== checking requires ===\n")
        # file.write("* size of d:" + str(floor(int(log(d)) / log(2))) + "\n")
        # file.write("* y0:" + str(y0) + "\n")
        # file.write("* Y:" + str(Y) + "\n")
        # file.write("* x0:" + str(x0) + "\n")
        # file.write("* X:" + str(X) + "\n")
        # file.write("* d < N**0.292?:" + str(d < N ** (0.292)) + "\n")
        # file.write("* |y| < Y?" + str(abs(y0) < Y) + "\n")
        # file.write("* |x| < X?" + str(abs(x0) < X) + "\n")

    # console
    if debug:
        print("=== running algorithm ===")
        start_time = time.time()

    solx, soly = boneh_durfee(pol, e, m, t, X, Y, file)

    if solx > 0:
        print("=== solution found ===")
        if True:
            print("x:", solx)
            print("y:", soly)

        d = int(pol(solx, soly) / e)
        print("private key found:", d)

        msg = 1145
        cipher = pow(msg, e, N)
        if msg != pow(cipher, d, N):
            print("private key is NOT CORRECT!")
        else:
            print("\n=== congratulations! we hacking it!===")

    print("=== %s seconds ===" % (time.time() - start_time))


def crt_rsa_example():
    # Configs
    file_addr = "debug_crt_rsa_bf.txt"
    file = open(file_addr, "w+")
    length_N = Integer(2048)
    length_d = RealNumber(0.28)

    # while length_d <= RealNumber(0.24):
        # for i in range(5, 5):  # size of the lattice (bigger, better,slower)
    delta = RealNumber(0.24)
    m = Integer(5)
    for i in range(0, 5):  # size of the lattice (bigger, better,slower)
        m+=Integer(1)
        print("m:",m)
        delta = RealNumber(0.24)
        while delta <= RealNumber(0.30):
            print("delta:",delta)
            crt_rsa_logic(file, length_N, length_d, delta, m)
            delta += RealNumber(0.01)
            print()
    # crt_rsa_logic(file, length_N, length_d, delta, m)
                # delta += RealNumber(0.01)
            # length_d += RealNumber(0.01)
    length_N = Integer(2048)
    length_d = RealNumber(0.27)
    m = Integer(5)
    delta = RealNumber(0.24)
    # crt_rsa_logic(file, length_N, length_d, delta, m)


def multiprime_logic(file, length_N, length_d, delta, m):
    print("m:", m)
    print("length_d:", length_d)
    p = Integer(next_prime(2 ** int(round(length_N / 3))))
    # q = Integer(next_prime(round(pi.n() * p)))
    q = Integer(next_prime(next_prime(p)))
    q = Integer(q)
    q2 = Integer(next_prime(q))
    N = Integer(p * q * q2)
    phiN = (p - 1) * (q - 1) * (q2 - 1)

    d = Integer(int(pow(N, length_d)))
    if d % 2 == 0:
        d += 1
    while gcd(d, phiN) != 1:
        d += 2
    e = d.inverse_mod(phiN)

    t = floor((1 - 2 * delta) * m)  # from Herrmann and May
    Y = floor(pow(N, RealNumber(2 / 3)))
    X = floor(2 * floor(N**delta)) 
    A = Integer((N + 1) / 2)

    # s & k in k(A+s) = 1 mod(e)
    y0 = floor(phiN - N - 1) / 4  # s
    x0 = floor((1 - e * d) / (A + y0)) + 1  # k

    # print("x0:",x0)

    # Problem put in equation
    P = PolynomialRing(
        ZZ,
        names=(
            "x",
            "y",
        ),
    )
    (
        x,
        y,
    ) = P._first_ngens(2)
    pol = x * (A + y) + 1

    file.write("p: 0x" + p.hex() + "\n\n")
    file.write("q: 0x" + q.hex() + "\n\n")
    file.write("r: 0x" + q2.hex() + "\n\n")
    file.write("d: 0x" + d.hex() + "\n\n")
    file.write("m: 0x" + m.hex() + "\n\n")
    file.write("t: 0x" + t.hex() + "\n\n")

    file.write("Result of " + str(length_d) + " and " + str(m) + "\n")
    file.write("* delta:" + str(delta) + "\n")

    file.write("* size of d:" + str(floor(int(log(d)) / log(2))) + "\n")
    file.write("* |y| < Y:" + str(abs(y0) < Y) + "\n")
    file.write("* |x| < X:" + str(abs(x0) < X) + "\n")
    file.write("* d < N**0.292?:" + str(d < N ** (0.292)) + "\n")

    print("=== running algorithm ===")
    start_time = time.time()

    solx, soly = boneh_durfee(pol, e, m, t, X, Y, file)

    if solx > 0:
        print("=== solution found ===")
        if True:
            print("x:", solx)
            print("y:", soly)

        d = int(pol(solx, soly) / e)
        print("private key found:", d)

        msg = 1145
        cipher = pow(msg, e, N)
        if msg != pow(cipher, d, N):
            print("private key is NOT CORRECT!")
            file.write("private key is NOT CORRECT!\n")
        else:
            file.write("=== congratulations! we hacking it!===\n")

    print("=== %s seconds ===" % (time.time() - start_time))
    file.write("=== %s seconds ===" % (time.time() - start_time) + "\n\n")


def multiprime_example():
    # Configs
    file_addr = "debug_Multiprime_bf.txt"
    file = open(file_addr, "w+")
    length_N = Integer(1000)
    length_d = RealNumber(0.17)

    m = Integer(7)
    delta = RealNumber(0.01)
    multiprime_logic(file, length_N, length_d, delta, m)
    
    return 
    
    while length_d <= RealNumber(0.17):
        for i in range(6, 9):  # size of the lattice (bigger, better,slower)
            m = Integer(i)
            print("m:",m)
            delta = RealNumber(0.13)
            while delta <= RealNumber(0.19):
                print("delta:",delta)
                multiprime_logic(file, length_N, length_d, delta, m)
                delta += RealNumber(0.01)
        length_d += RealNumber(0.01)


if __name__ == "__main__":
    debug = True
    multiprime_example()
    # crt_rsa_example()
