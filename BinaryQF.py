from typing import Any as Any, Generator, NoReturn as NoReturn
from sage.calculus.var import var
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.arith.misc import is_square, is_prime, next_prime, gcd
from sage.misc.functional import isqrt, sqrt
from sage.functions.other import floor, ceil
from sage.matrix.constructor import Matrix


# The constructor
#
# binaryQF([a, b, c])
#    
# This is the constructor of the class 'binaryQF'. 
# a, b, c are the integer coefficients of the quadratic form
# q(x, y) = ax^2 + bxy + cy^2. 
#
# The main function
#
# Return all numbers represented by the quadratic form subject to the
# constraint _subset_, which is _all_, _primitively_ or _prime_ up to 
# the bound _upto_. 
#
# represented_positives( 
#    upto,    # search the range (1..upto)
#    subset,  # 'all' or 'primitively' or 'prime', default 'all'
#    verbose  # print messages, default True
# )
    

class binaryQF():
    """
    A binary quadratic form over Z.
    Input: a list of 3 entries: [a, b, c]
    """

    def __init__(self, 
        abc: list[int]
    ) -> None:

        self._a: Integer = ZZ(abc[0])
        self._b: Integer = ZZ(abc[1])
        self._c: Integer = ZZ(abc[2])


    def discriminant(self
    ) -> Integer:

        return self._b ** 2 - 4 * self._a * self._c


    def sqr_disc(self, 
        M: int, 
        primitively: bool = False
    ) -> list[Integer]:

        d = self.discriminant()
        if d == 0:
            raise ValueError("Discriminant must not be zero")

        a, b, c = self._a, self._b, self._c
        if a == 0 and c == 0: # and b != 0
            return [b * n for n in range(1, 1 + M // abs(b))]

        D = Integer(d).sqrtrem()[0]

        # a must be != 0
        if a == 0: # then c <> 0; swap
            a = c
            c = 0

        k = 2 * D; m = 4 * a * D
        u = b + D; v = b - D
        S = set[Integer]()

        # Solvability in Z.
        for n in range(1, M + 1):
            h = 4 * a * n  # a <> 0 and n <> 0
            for t in h.divisors():
                g = h // t
                if k.divides(g - t) and m.divides(g * u - t * v):

                    if primitively:
                        y = (g - t) // k
                        x = var('x')
                        eq = a * x * x + b * x * y + c * y * y
                        R = (eq - n).roots(multiplicities = False, ring = ZZ)
                        x = R[0]

                        if gcd(x, y) == 1:
                            S.add(n)
                            break
                    else:
                        S.add(n)
                        break

        return sorted(list(S))


    def imag_prime(self, 
        M: int
    ) -> list[Integer]:

        solve = pari('qfbsolve')
        Q = pari('Qfb')(self._a, self._b, self._c)
        p = 1
        r = []

        while True:
            p = next_prime(p)
            if p > M: break
            if solve(Q, p):
                r.append(p)
        return r


    def imag_primitively(self, 
        M: int
    ) -> list[Integer]:

        a, b, c = self._a, self._b, self._c
        d = c - b * b / (4 * a) 
        A : list[Integer] = []

        for y in range(1 + isqrt(M / d)):
            r = y * b / (2 * a)
            s = sqrt((M - d * y * y) / a)
            for x in range(ceil(-s -r), 1 + floor(s - r)):
                if gcd(x, y) == 1:
                    A.append(a * x^2 + b * x * y + c * y^2) 

        return sorted(list(set(A)))


    def imag_all(self, 
        M: int
    ) -> list[Integer]:

        L = [2 * ZZ(self._a), ZZ(self._b), ZZ(self._b), 2*ZZ(self._c)]
        G = Matrix(ZZ, 2, 2, L)
        A = pari('qfrep')(G, M, 1)
        return [k + 1 for k in range(M) if A[k] > 0]


    def _primitive_reps(self, 
        a: Integer, 
        h: Integer, 
        b: Integer, 
        M: int, 
        S: set[Integer]
    ) -> None:

        if a <= M :
            S.add(a)
            if b <= M :
                S.add(b)
                if a <= (M - b) and h <= (M - a - b) :
                    if a <= (M - a - h) :
                        self._primitive_reps(a, h + 2 * a, a + b + h, M, S)
                    if b <= (M - b - h) :
                        self._primitive_reps(a + b + h, h + 2 * b, b, M, S)


    def positive_primitives(self, 
        M: int, 
        primitively: bool
    ) -> set[Integer]:

        a, b, c = self._a, self._b, self._c

        S = set[Integer]()
        while True:
            new_val = a + b + c
            if new_val > 0 :
                self._primitive_reps(a, b + 2 * a, new_val, M, S)
                b += 2 * c
                a = new_val
            elif new_val < 0 :
                b += 2 * a
                c = new_val
            if a == self._a and b == self._b and c == self._c:
                break

        if not primitively :
            X = set[Integer]()
            for p in S:
                q = t = 1
                while q <= M :
                    X.add(q)
                    q = t * t * p
                    t += 1
            S = X

        return S


    def reduce_real(self
    ) -> list[Integer]:

        d = self.discriminant()
        if is_square(d):
            raise ValueError("Form must not have square discriminant")

        droot = Integer(d).sqrtrem()[0]
        a, b, c = self._a, self._b, self._c

        while a <= 0 or c >= 0 or b <= abs(a + c):

            cAbs = c
            if cAbs < 0: cAbs *= -1

            # cAbs = 0 will not happen for a non square form
            delta = (b + droot) // (2 * cAbs)
            if c < 0: delta *= -1
            aa = c
            bb = 2 * c * delta - b
            cc = c * delta * delta - b * delta + a
            a, b, c = aa, bb, cc

        return [a, b, c]


    def reduce_imag(self
    ) -> list[Integer]:

        a, b, c = self._a, self._b, self._c
        if a < 0: a, b, c = -a, -b, -c
        d = self.discriminant()

        while True:
            A = ( a == c and b < 0) or (c < a)
            B = (-a == b and a < c) or (a < abs(b))

            if not (A or B) : break

            if A: a, b, c = c, -b, a

            if B:
                b -= 2 * a * (b // (2 * a))
                if abs(b) > a: b -= 2 * a
                c = (b * b - d) // (4 * a)

        return [a, b, c]


    def is_reduced(self
    ) -> bool:

        a, b, c = self._a, self._b, self._c
        return (-a < b <= a < c) or (ZZ(0) <= b <= a == c)


    def reduced_form(self):
        """
        Returns the unique reduced form equivalent to binaryQF(a, b, c)
        """

        if self.is_reduced() :
            return self

        if self.discriminant() >= 0:
            r = self.reduce_real()
        else:
            r = self.reduce_imag()

        return binaryQF(r)


    def represented_positives(self, 
        upto: int, 
        subset: str = "all", 
        verbose: bool = True
    ) -> list[Integer]:
        """
        subset = "all" or "primitively" or "prime"
        """

        prime = False or subset == "prime"
        primitively = False or subset == "primitively"

        d = self.discriminant()
        if d == 0:
            raise ValueError("discriminant must not be 0")

        a, b, c = self._a, self._b, self._c
        if verbose:
            print("Original form ", [a, b, c], "with discriminant", d)

        if is_square(d):
            if verbose:
                print("Square discriminant!")

            if prime:  # for efficiency
                primitively = False
            pp = self.sqr_disc(upto, primitively)
            if prime:
                pp = list(filter(is_prime, pp))

        else:

            R = self.reduced_form()
            if verbose:
                print("Reduced form  ", [R._a, R._b, R._c])

            if d < 0:

                if prime:
                    pp = R.imag_prime(upto)
                else:
                    if primitively:
                        pp = R.imag_primitively(upto)
                    else:
                        pp = R.imag_all(upto)

            # real case, indefinite form
            else: # d > 0 and not square

                if prime:  # for efficiency
                    primitively = True
                pp = R.positive_primitives(upto, primitively)
                if prime:
                    pp = list(filter(is_prime, pp))
                pp = sorted(pp)

        if verbose:
            msg0 = "primes" if prime else "positive integers"
            msg1 = "primitively" if primitively else ""
            msg2 = "represented up to"
            print("There are", len(pp), msg0, msg1, msg2, upto)

        return pp

# The OEIS-query function</h2>
#  
# oeis_bqf(q, filter, upto, terse)
# 
# The function tries to find sequences in the OEIS whose terms are represented
# by the binary quadratic form with coefficients $q = [a, b, c]$ and which are
# restricted according to the _filter_, which is one of _all_, _primitively_,
# _prime_ or _tutti_. The parameter _upto_ gives the upper bound of the search 
# range, which is 100 by default. If _terse_ is _True_ the output will be a 
# one-liner; otherwise the output is more verbose. With the parameter _values_ 
# you can switch off the display of the values; it is set to _True_ by default. 
# To use the function you have to be connected to the Internet.

def oeis_bqf(
    abc: list[int], 
    upto: int = 100,
    filter: str = 'all',
    terse: bool = True,
    values: bool = True
):
    
    if filter == 'tutti':
        oeis_bqf([1, 1, 1], upto, 'all', terse=False) 
        oeis_bqf([1, 1, 1], upto, 'primitively', terse=False) 
        oeis_bqf([1, 1, 1], upto, 'prime', terse=False) 
        return
        
    reps = []
    Q = binaryQF(abc)
    reps = Q.represented_positives(upto, filter, verbose = not terse)

    d = abc[1] ** 2 - 4 * abc[0] * abc[2]
    if reps == []:
        print(f"No the representatives below {upto}.")
        print([d], abc, filter)
        return

    reps = reps[:min(20, upto)]
    if values and not terse: print(reps) 
    search = oeis(reps, 4)

    found = []
    if search != []:
        if not terse: print(search)
        found = [seq.id() for seq in search] 
        if not terse: 
            if found == []: print("No sequence found in the OEIS.")
            else: print(found)
        if terse: 
            if found == []: 
                print([d], abc, filter)
                print("No sequence found in the OEIS.")
            else: 
                print([d], abc, filter, found)
            if values: 
                print(reps) 
                print()
        else: print()
