def modpow(base, exponent, modulus):
    if modulus == 1:
        return 0
    # Assert :: (modulus - 1) * (modulus - 1) does not overflow base
    result = 1
    base = base % modulus
    while exponent > 0:
        if exponent & 1 == 1:
           result = (result * base) % modulus
        exponent >>= 1
        base = (base * base) % modulus
    return result

def modinv(x, p):
    """
    Returns the inverse in the multiplicative group of Z_p.
    If none exists, which occurs if the gcd(x, p) != 1, then None is returned.
    """
    result = pow(x, p - 2, p)
    # p was not a prime number -> inverse may not have been found
    return None if x * result % p != 1 else result 

def mulinv(x, mod = None):
    """
    Returns modinv(x, mod) if a mod is given, else 1 / x.
    """
    return 1 / x if mod == None else modinv(x, mod)
    
def eulerphi(n, primes=None):
    # calculates Euler's totient function for a given number n in O(n)
    # optionally, n's prime factors may be given which reduces the complexity to O(len(primes))
    if primes is None:
        primes = primefacs(n).keys()
    return n * product(p - 1 for p in primes) // product(p for p in primes) 
