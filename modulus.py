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
