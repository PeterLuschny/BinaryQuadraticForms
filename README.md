# Binary Quadratic Forms ...
 
Positive Numbers represented by a binary quadratic form.

This SageMath notebook implements methods for calculating the numbers represented by a binary quadratic form in a uniform and efficient way and hides the complexity of the situation by providing a simple interface for the user.


The original SageMath notebook (Python2, 2014) was documented at OEIS: [BinaryQuadraticForms-OEIS](https://oeis.org/wiki/User:Peter_Luschny/BinaryQuadraticForms)

The updated notebook (Python 3.10, 2024) is at GitHub, as well as the pure Python code from the updated notebook: [BinaryQuadraticForms-GitHub](https://github.com/PeterLuschny/BinaryQuadraticForms)

A HTML version of the notebook can be found at: [BinaryQuadraticForms-HTML](https://luschny.de/math/seq/binaryqf/BinaryQF.html)


# ... at your fingertips:

Go to [SageMath CellServer](https://sagecell.sagemath.org/). Enter one example after the other:

load('https://raw.githubusercontent.com/PeterLuschny/BinaryQuadraticForms/main/BinaryQF.sage')


// --- Example 1:

Q = binaryQF([1, 0, -2])

Q.represented_positives(30)


// --- Example 2:

oeis_bqf([1, 1, -1], 100, 'primitively', terse=False) 


// --- Example 3:

oeis_bqf([1, 1, 1], 100, 'tutti')
