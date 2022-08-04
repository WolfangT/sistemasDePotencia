#!/usr/bin/env python3
# comun.py
# Herramientas de Wolfang para calculos rapidos en sistemas de potencia


# Standard Library
import cmath
from cmath import pi

# Pip
import numpy as np

# Conf librerias

np.set_printoptions(precision=2, suppress=True, linewidth=180)

# Constantes

G = 10 ** 9
M = 10 ** 6
K = 10 ** 3
m = 10 ** -3
u = 10 ** -6
n = 10 ** -9

inf = float("inf")

Vb13 = 13800
Vb24 = 24000


# Funciones de ayuda


def paralelo(elementos):
    """Calcula el inverso de la suma de los inversos de una serie de elementos"""
    return 1 / sum(1 / e for e in elementos)


def degrees(rad):
    """converts radians to degress"""
    return rad * 180 / pi


def radians(degs):
    """converts degrees to radians"""
    return degs * pi / 180


def rect(mag, deg):
    """creates a complex number from magnitude and degrees"""
    return cmath.rect(mag, radians(deg))


def ang(val):
    """regresa el angulo en grados"""
    mag, rad = cmath.polar(val)
    return degrees(rad)


def display_single(val):
    return f"{round(abs(val), 6):11}"


def display_polar(cmplx, unit=""):
    mag, rad = cmath.polar(cmplx)
    deg = degrees(rad)
    return f"{mag:11.6f} ∠ {deg:7.2f}° {unit}"


def display_rect(cmplx, u_real="", u_imag="J", u_fin=""):
    return f"{cmplx.real:11.6f} {u_real} + {cmplx.imag:11.6f} {u_imag} {u_fin}"


def cambio_base(val, Sn, Vn, Sb, Vb):
    """Cabio de Base de valores nominales a valores bases"""
    return val * ((Vn / Vb) ** 2) * (Sb / Sn)


# Clases de Valores


class Valor(complex):
    TIPOS = {
        "S": (("W", "VAR", "VA"), "rect"),
        "V": ((None, None, "V"), "polar"),
        "I": ((None, None, "A"), "polar"),
        "Z": ((None, None, "Ω"), "rect"),
        "Y": ((None, None, "℧"), "rect"),
    }
    TIPO = None

    def __new__(cls, name, comp, barra=None, muestra=None):
        return super().__new__(cls, comp)

    def __init__(self, name, comp, barra=None,  muestra=None):
        self.name = name
        self.unidad = self.TIPOS[self.TIPO][0] if self.TIPO else ("", " J", "")
        self.muestra = muestra or self.TIPOS[self.TIPO][1] if self.TIPO else "rect"
        self.barra = barra
        super().__init__()

    def __str__(self):
        u_fin = self.unidad[2] + (" pu" if self.barra else "")

        if self.muestra == "rect":
            if self.unidad[:2] == (None, None):
                val = display_rect(self, u_fin=u_fin)
            else:
                u_r = self.unidad[0] + (" pu" if self.barra else "")
                u_i = self.unidad[1] + (" pu" if self.barra else "")
                val = display_rect(self, u_r, u_i)
        elif self.muestra == "polar":
            val = display_polar(self, u_fin)
        return f"{self.name}: {val}"

    def en_real(self):
        if self.barra:
            base = getattr(self.barra, self.TIPO + "b")
            valor = self * base
        else:
            valor = self
        return self.__class__(self.name, valor).__str__()


class S(Valor):
    TIPO = "S"


class V(Valor):
    TIPO = "V"


class I(Valor):  # noqa: E742
    TIPO = "I"


class Z(Valor):
    TIPO = "Z"


class Y(Valor):
    TIPO = "Y"
