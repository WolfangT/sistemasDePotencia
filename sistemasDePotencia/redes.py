#!/usr/bin/env python3
# redes.py
# Herramientas de Wolfang para calculos rapidos en sistemas de potencia 1


import numpy as np
from sistemasDePotencia.potencia import Elemento, Barra, GeneradorIdeal, CargaIdeal
from sistemasDePotencia.despacho import GeneradorDespacho, CargaDespacho
from sistemasDePotencia.comun import inf, display_polar, rect, ang


# Red de Potencia


class RedPotencia:
    """Obtiene una matris de admitancias"""

    def __init__(self, barras: list[Barra], elementos: list[Elemento], generadores: list[GeneradorIdeal]):
        self.barras = barras
        self.elementos = elementos
        self.generadores = generadores
        self.calcularMatris()

    @property
    def matrisAdmitancias(self):
        return self._matris

    @property
    def matrisImpedancias(self):
        return np.linalg.inv(self._matris)

    @property
    def barrasGeneracion(self):
        return self._barras_gen

    @property
    def barrasCarga(self):
        return self._barras_carga

    @property
    def barrasTransmision(self):
        return self._barras_trans

    def calcularMatris(self, generadores: list[str], cargas: list[CargaIdeal] = []):
        nb = len(self.barras)
        self._matris = np.zeros((nb, nb), dtype=np.complex64)
        self._barras_gen = set()
        self._barras_carga = set()
        self._barras_trans = set()

        for elemento in self.elementos + generadores + cargas:
            bp = self.barras.index(elemento.bp)
            try:
                bs = self.barras.index(elemento.bs)
            except AttributeError:
                bs = None
            if bs is None:  # elemento entre barra y tierra
                self._matris[bp][bp] += elemento.Y
                if isinstance(elemento, CargaIdeal):
                    self._barras_carga.add(elemento.bp)
                if isinstance(elemento, GeneradorIdeal):
                    self._barras_gen.add(elemento.bp)
            else:  # elemento entre dos barras
                self._matris[bp][bp] += elemento.Yp0 + elemento.Yps
                self._matris[bs][bs] += elemento.Ys0 + elemento.Yps
                self._matris[bp][bs] -= elemento.Yps
                self._matris[bs][bp] -= elemento.Yps
        for barra in self.barras:
            if barra not in self._barras_carga and barra not in self._barras_gen:
                self._barras_trans.add(barra)


class GaussSiedel:
    """Ejecuta Gauss Siedel para resolver un sistema de potencia

    Requiere una lista de barras, Potencias activas, reactivas,
    un factor de aceleracion y la matris de admitancias.

    La barra 1 siempre es la de compensacion
    """

    TIPOS = {1: "Compensación", 2: "Volt. Regulado", 3: "P - Q"}

    def __init__(self, barras, Vn, Pn, Qn, alpha, Ym):
        self.barras = barras
        self.Vn = Vn
        self.Qn = Qn
        self.Pn = Pn
        self.tipos = [1] + [2 if v else 3 for v in Vn[1:]]
        self.alpha = alpha
        self.Ym = Ym
        self.iter = 0
        self.V = [[v or 1 + 0j for v in self.Vn]]
        self.P = [p or 0 for p in self.Pn]
        self.Q = [[q or 0j for q in self.Qn]]
        self.error = [[inf for v in self.Vn]]

    def __str__(self):
        k = self.iter
        lineas = "\n".join(
            [
                " | ".join(
                    [
                        f"{barra.nombre:3}",
                        f"{display_polar(v)}",
                        f"{p.real:11.06f}" if p else "-".center(11),
                        f"{q.imag:11.06f}" if q else "-".center(11),
                        f"{self.TIPOS[tipo]}",
                    ]
                )
                for (barra, v, p, q, tipo) in zip(self.barras, self.V[k], self.P, self.Q[k], self.tipos)
            ]
        )
        return (
            f"Metodo de Gauss Siedel - iter {self.iter}: (α={self.alpha}, error={self.max_error*100:4.2f}%)\n"
            + " | ".join([" # ", "V".center(23), "P".center(11), "Q".center(11), " Tipo "])
            + f"\n{lineas}\n"
        )

    @property
    def max_error(self):
        return max(self.error[-1])

    def calcQ(self, k, i):
        return (
            -(
                self.V[k - 1][i].conjugate()
                * (
                    sum(self.Ym[i][j] * self.V[k][j] for j in range(i))
                    + sum(self.Ym[i][j] * self.V[k - 1][j] for j in range(i, len(self.barras)))
                )
            ).imag
            * 1j
        )

    def calcV(self, k, i, Q=None):
        Q = Q or self.Qn
        return (
            1
            / self.Ym[i][i]
            * (
                ((self.P[i] - Q[i]) / self.V[k - 1][i].conjugate())
                - sum(self.Ym[i][j] * self.V[k][j] for j in range(i))
                - sum(self.Ym[i][j] * self.V[k - 1][j] for j in range(i + 1, len(self.barras)))
            )
        )

    def acel(self, Vv, Vn):
        return Vv + self.alpha * (Vn - Vv)

    def iteracion(self):
        self.iter += 1
        k = self.iter
        self.V.append([self.V[k - 1][i] if self.tipos[i] in (1, 2) else None for i, b in enumerate(self.barras)])
        self.Q.append([self.Q[k - 1][i] if self.tipos[i] == 3 else None for i, b in enumerate(self.barras)])
        self.error.append([0 for b in self.barras])
        for i, barra in enumerate(self.barras):
            if self.tipos[i] == 1:
                continue
            elif self.tipos[i] == 2:
                self.Q[k][i] = self.calcQ(k, i)
                self.Q[k][i] = self.acel(self.Q[k - 1][i], self.Q[k][i])
                self.V[k][i] = rect(abs(self.V[k][i]), ang(self.calcV(k, i, self.Q[k])))
            elif self.tipos[i] == 3:
                self.V[k][i] = self.calcV(k, i)
                self.V[k][i] = self.acel(self.V[k - 1][i], self.V[k][i])
            self.error[k][i] = abs((self.V[k][i] - self.V[k - 1][i]) / self.V[k - 1][i])
        return str(self)

    def resolver(self, error=0.01, max_iter=100):
        """Resolver el sistemas hasta obtener un error maximo"""
        for i in range(max_iter):
            self.iteracion()
            res = [e < error for e in self.error[-1]]
            if all(res):
                break
        return str(self)

    @property
    def matrisVoltajes(self):
        return np.array([[v] for v in self.V[-1]])


def reduccionKron(Ym, p):
    """Elimina un nodo p de la matris"""
    p -= 1
    x, y = Ym.shape
    if x != y:
        raise ValueError("La matris debe ser cuadrada")
    if not 0 <= p <= x:
        raise ValueError("El nodo a eliminar debe estar dentro de la matris")
    Yn = np.array([[0 + 0j] * (x - 1)] * (x - 1))
    a = 0
    for i in range(x):
        if i == p:
            continue
        b = 0
        for j in range(y):
            if j == p:
                continue
            Yn[a][b] = Ym[i][j] - (Ym[i][p] * Ym[p][j] / Ym[p][p])
            b += 1
        a += 1
    return Yn


class RedSimplificada:
    """Red simplificada de barras de despacho

    Determina las matrices de transformacion C
    """

    def __init__(self, *barras):
        self.barras = barras
        self.gens = [gen for barra in self.barras for gen in barra.gens]
        self.cargas = [barra.carga for barra in self.barras if barra.carga]
        self._C12 = None
        self._C34 = None

    @property
    def corrientes1(self):
        return np.array([[f"I{barra.nombre}"] for barra in self.barras])

    @property
    def corrientes2(self):
        return np.array([[f"I{gen.nombre}"] for gen in self.gens] + [[f"I{carga.nombre}"] for carga in self.cargas])

    @property
    def corrientes3(self):
        return np.array([[f"I{gen.nombre}"] for gen in self.gens] + [["Ist"]])

    @property
    def corrientes4(self):
        return np.array([[f"I{gen.nombre}"] for gen in self.gens] + [["IT"]])

    @property
    def potenciasG(self):
        return np.array([[f"P{gen.nombre}"] for gen in self.gens] + [["1"]])

    @property
    def C12(self):
        if self._C12 is None:
            self._C12 = self.separacionDeBarras()
        return self._C12

    @property
    def C34(self):
        if self._C34 is None:
            self._C34 = self.cambioDeReferencia()
        return self._C34

    def separacionDeBarras(self):
        """Separa las barras en barras de cargas y de generación"""
        cols = self.gens + self.cargas
        C12 = np.zeros((len(self.barras), len(cols)), dtype=np.complex64)
        for y, barra in enumerate(self.barras):
            for x, col in enumerate(cols):
                if isinstance(col, GeneradorDespacho):
                    if col in barra.gens:
                        C12[y][x] = 1
                elif isinstance(col, CargaDespacho):
                    if col is barra.carga:
                        C12[y][x] = 1
        return C12

    def unificacionDeCargas(self, Is):
        """Unifica las cargas segun las corrientes que consumen

        Args:
            Is (list[float]): lista de la corrientes que consumen cada carga
        """
        ng = len(self.gens)
        nc = len(self.cargas)
        if len(Is) != len(self.cargas):
            raise ValueError("Se requiere de una corriente por carga")
        It = sum(Is)
        k = [i / It for i in Is]
        C23 = np.zeros((ng + nc, ng + 1), dtype=np.complex64)
        for y in range(ng):
            C23[y][y] = 1
        for y, i in enumerate(k):
            C23[ng + y][ng] = i
        return C23

    def cambioDeReferencia(self):
        """Cambio de referencia de carga a tierra"""
        ng = len(self.gens)
        C34 = np.identity(ng + 1, dtype=np.complex64)
        for x in range(ng):
            C34[ng][x] = -1
        return C34

    def cambioDeVariable(self, IT, ms):
        """Cambio de variable de corriente a potencia

        Args:
            IT (float): Corriente de tierra
            ms (list[float]): lista de factores de  trasnformacion
        """
        ng = len(self.gens)
        C4G = np.identity(ng + 1, dtype=np.complex64)
        for n, m in enumerate(ms):
            C4G[n][n] = m
        C4G[ng][ng] = IT
        return C4G
