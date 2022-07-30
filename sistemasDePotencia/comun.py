#!/usr/bin/env python3
# comun.py
# Herramientas de Wolfang para calculos rapidos en sistemas de potencia


# Standard Library
import cmath
from cmath import pi

# Constantes

M = 10**6
K = 10**3
m = 10**-3
u = 10**-6

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


def radians(degrees):
    """converts degrees to radians"""
    return degrees * pi / 180


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
