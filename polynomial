#!/usr/bin/python3

import math

import modulus

from modulus import *

def longer(a, b):
    return a if len(a) >= len(b) else b

def shorter(a, b):
    return a if len(a) < len(b) else b

class Polynomial:
    
    def degree(self):
        return len(self.cof) - 1
    
    def copy(self):
        return Polynomial(self.cof.copy())
    
    def trim(self):
        for i in range(len(self.cof)):
            if self.cof[i] != 0:
                return Polynomial(self.cof[i:])
        # all coefficients are zero
        return Polynomial([0])
    
    def cmod(self, field_size):
        """
        Applies a modulo operation to all coefficients.
        """
        result = []
        for i in range(len(self.cof)):
            result.append(self.cof[i] % field_size)
        return Polynomial(result)
    
    def pretty_repr(self):
        d = self.degree()
        result = []
        for i in range(len(self.cof)):
            coeff = self.cof[i]
            if coeff != 0:
                result.append((str(coeff) if coeff != 1 or d == 0 else "") + ("x^" + str(d) if d > 1 else "x" if d == 1 else ""))
            d -= 1
            
        result_str = result[0]
        for i in range(1, len(result)):
            ele = result[i]
            result_str += (" - " + ele[1:]) if ele[0] == '-' else (" + " + ele)
        return result_str
    
    def __init__(self, *args):
        if len(args) == 1:
            args = args[0]
            if isinstance(args, int) or isinstance(args, float):
                self.cof = [args]
            elif not isinstance(args, list):
                raise AssertionError("coefficients must be a list")
        elif not args:
            raise AssertionError("coefficient list can not be empty")
        self.cof = list(args)
        
    def __str__(self):
        return "P" + str(self.cof)
    
    def __repr__(self):
        return str(self)
        
    def __getitem__(self, degree):
        """
        Returns the coefficient of the variable with the given exponent.
        """
        if degree < 0:
            raise IndexError("negative degrees not allowed")
        if degree > self.degree():
            raise IndexError("degree too large")
        return self.cof[self.degree() - degree]
    
    def __setitem__(self, degree, value):
        """
        Sets the coefficient of the variable with the given exponent.
        """
        if degree < 0:
            raise IndexError("negative degrees not allowed")
        deg_diff = degree - self.degree()
        if deg_diff > 0:
            self.cof = [0] * deg_diff + self.cof
            self.cof[0] = value
        else:
            self.cof[self.degree() - degree] = value
        return self
    
    def __neg__(self):
        result = [0] * len(self.cof)
        for i in range(len(self.cof)):
            result[i] = -self.cof[i]
        return Polynomial(result)
        
    def __add__(self, pol):
        if isinstance(pol, int) or isinstance(pol, float):
            result = self.copy()
            result[0] += pol
            return result
        if not isinstance(pol, Polynomial):
            raise AssertionError("can only add with Polynomial")
        lng = longer(self.cof, pol.cof)
        srt = shorter(self.cof, pol.cof)
        deg = len(lng)
        result = []
        diff = len(lng) - len(srt)
        
        for i in range(deg):
            i_s = i - diff
            coeff = lng[i] if i_s < 0 else lng[i] + srt[i_s]
            #print(str(coeff) + " = " + str(lng[i]) + " + " + str(srt[i_s]))
            result.append(coeff)
        return Polynomial(result)
    
    def __sub__(self, pol):
        return self + (-pol)
    
    def __mul__(self, pol):
        if isinstance(pol, int) or isinstance(pol, float):
            result = self.copy()
            for i in range(len(self.cof)):
                result[i] *= pol
            return result
        if not isinstance(pol, Polynomial):
            raise AssertionError("can only multiply with Polynomial")
        deg = self.degree() + pol.degree()
        result = Polynomial([0] * (deg + 1))
        
        for i in range(len(self.cof)):
            for j in range(len(pol.cof)):
                result[i + j] += self[i] * pol[j]
        return result
    
    def divmod(self, pol, gf = None):
        if not isinstance(pol, Polynomial):
            raise AssertionError("can only divide by Polynomial")
        quo = Polynomial([0])
        tmp = self.trim()
        pol = pol.trim()
        while tmp.degree() >= pol.degree():
            #print("tmp = " + str(tmp) + " with degree " + str(tmp.degree()))
            coeff = tmp[tmp.degree()] * mulinv(pol[pol.degree()], gf)
            if gf != None:
                coeff %= gf
            deg = tmp.degree() - pol.degree()
            
            quo[deg] = coeff
            #print("quo = " + str(quo))
            tmpFac = Polynomial([0])
            tmpFac[deg] = coeff
            tmp -= pol * tmpFac
            #print("product = " + str(pol * tmpFac))
            #print("new_tmp = " + str(tmp))
            # zero out to prevent loop caused by float imprecisions
            tmp[tmp.degree()] = 0
            if gf != None:
                tmp = tmp.cmod(gf)
            tmp = tmp.trim()
        return [quo.trim(), tmp]

    def add(self, pol, gf = None):
        result = self + pol
        return result if gf == None else result.cmod(gf)
    
    def sub(self, pol, gf = None):
        result = self - pol
        return result if gf == None else result.cmod(gf)
    
    def mul(self, pol, gf = None):
        result = self * pol
        return result if gf == None else result.cmod(gf)
    
    def div(self, pol, gf = None):
        return self.divmod(pol, gf)[0]
    
    def mod(self, pol, gf = None):
        return self.divmod(pol, gf)[1]
    
    def modpow(self, exponent, modulus, gf = None):
        """
        Exponentiates this polynomial.
        """
        if exponent < 0:
            raise ArithmeticError("negative exponents not supported")
        result = Polynomial([1])
        for i in range(exponent):
            result = result.mul(self, gf).mod(modulus, gf)
        return result
    
    def __truediv__(self, pol):
        return self.divmod(pol)[0]
    
    def __mod__(self, pol):
        return self.divmod(pol)[1]
    
    def __pow__(self, exponent):
        """
        Exponentiates this polynomial.
        """
        if exponent < 0:
            raise ArithmeticError("negative exponents not supported")
        result = Polynomial([1])
        for i in range(exponent):
            result *= self
        return result
    
    def __iadd__(self, pol):
        self.cof = (self + pol).cof
        return self
    
    def __isub__(self, pol):
        self.cof = (self - pol).cof
        return self
    
    def __imul__(self, pol):
        self.cof = (self * pol).cof
        return self
    
    def __imod__(self, pol):
        self.cof = (self % pol).cof
        return self
    
    def __ipow__(self, exponent):
        self.cof = (self**exponent).cof
    
if __name__ == "__main__":
    p = Polynomial([1, 2])
    print(p)
    print(-p)
    
    print(p + p)
    print(p - p)
    
    print("const tests")
    print(p + Polynomial([-1, -2]))
    print(p - Polynomial([1, 2]))
    p -= Polynomial([1, 2])
    print(p)
