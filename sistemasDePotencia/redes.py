#!/usr/bin/env python3
# redes.py
# Herramientas de Wolfang para calculos rapidos en sistemas de potencia 1

# Pip
import numpy as np
from texttable import Texttable

# Sistemas de Potencia
from sistemasDePotencia.comun import ang, inf, rect, display_polar
from sistemasDePotencia.despacho import Escenario, Combinacion, DespachoEconomico, GrupoCombinaciones
from sistemasDePotencia.potencia import Barra, Elemento, CargaIdeal, GeneradorIdeal

# Red de Potencia


class GaussSiedel:
    """Ejecuta Gauss Siedel para resolver un sistema de potencia

    Requiere una lista de barras, Potencias activas, reactivas,
    un factor de aceleracion y la matris de admitancias.

    La barra 1 siempre es la de compensacion
    """

    TIPOS = {1: "Compensación", 2: "Volt. Regulado", 3: "P - Q"}

    def __init__(self, barras, Vn, Pn, Qn, Ym, alpha=1):
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


class RedPotencia:
    """Obtiene una matris de admitancias"""

    def __init__(
        self,
        barras: list[Barra],
        elementos: list[Elemento],
        generadores: list[GeneradorIdeal],
        cargas: list[CargaIdeal],
    ):
        self.barras = barras
        self.elementos = elementos
        self.generadoresPosibles = generadores
        self.cargasPosibles = cargas
        self._generadores = None
        self._cargas = None
        self._matris = None
        self._barras_gen = None
        self._barras_carga = None
        self._barras_trans = None
        self._GS = None

    @property
    def matrisAdmitancias(self):
        return self._matris

    @property
    def barrasGeneracion(self):
        return self._barras_gen

    @property
    def barrasCarga(self):
        return self._barras_carga

    @property
    def barrasTransmision(self):
        """ "Barras sin cargas o generación"""
        self._barras_trans = []
        if self._barras_gen is not None and self._barras_carga is not None:
            for barra in self.barras:
                if barra not in self._barras_carga and barra not in self._barras_gen:
                    self._barras_trans.append(barra)
        return self._barras_trans[::-1]

    @property
    def generadores(self):
        """Los generadores que estaran activos segun la combinación"""
        return self._generadores

    @generadores.setter
    def generadores(self, combinacion):
        self._barras_gen = set()
        self._generadores = []
        for gen in self.generadoresPosibles:
            if combinacion is None or gen in combinacion:
                self._generadores.append(gen)
                self._barras_gen.add(gen.bp)

    @property
    def cargas(self):
        """Las cargas con sus valores segun la situación"""
        return self._cargas

    @property
    def GaussSiedel(self):
        return self._GS

    @property
    def Vg(self):
        if not self.GaussSiedel:
            return None
        return [self.GaussSiedel.V[-1][self.barras.index(gen.bp)] for gen in self.generadores]

    @property
    def Ig(self):
        if not self.GaussSiedel:
            return None
        return [(gen.S / v).conjugate() for gen, v in zip(self.generadores, self.Vg)]

    @property
    def Igt(self):
        if not self.GaussSiedel:
            return None
        return sum(self.Ig)

    @property
    def Vs(self):
        if not self.GaussSiedel:
            return None
        return [self.GaussSiedel.V[-1][self.barras.index(carga.bp)] for carga in self.cargas]

    @property
    def Is(self):
        if not self.GaussSiedel:
            return None
        return [(carga.S / v).conjugate() for carga, v in zip(self.cargas, self.Vs)]

    @property
    def Ist(self):
        if not self.GaussSiedel:
            return None
        return sum(self.Is)

    @property
    def It(self):
        if not self.GaussSiedel:
            return None
        return self.Igt + self.Ist

    @cargas.setter
    def cargas(self, escenario):
        self._barras_carga = set()
        self._cargas = []
        for carga in self.cargasPosibles:
            if escenario is not None and carga in escenario:
                carga.S = escenario[carga]
                self._cargas.append(carga)
                self._barras_gen.add(carga.bp)
            elif escenario is None:
                self._cargas.append(carga)
                self._barras_gen.add(carga.bp)

    def procesar(self, combinacion: Combinacion = None, escenario: Escenario = None):
        """Procesa todos los elementos, creando la matris de admitancia"""
        self._GS = None
        nb = len(self.barras)
        self._matris = np.zeros((nb, nb), dtype=np.complex64)
        self.generadores = combinacion
        self.cargas = escenario
        # procesar todos los elementos de la red
        for elemento in self.elementos + self.generadores + self.cargas:
            bp = self.barras.index(elemento.bp)
            try:
                bs = self.barras.index(elemento.bs)
            except AttributeError:
                bs = None
            if bs is None:  # elemento entre barra y tierra
                self._matris[bp][bp] += elemento.Y
            else:  # elemento entre dos barras
                self._matris[bp][bp] += elemento.Yp0 + elemento.Yps
                self._matris[bs][bs] += elemento.Ys0 + elemento.Yps
                self._matris[bp][bs] -= elemento.Yps
                self._matris[bs][bp] -= elemento.Yps

    def resolver(self, alpha=1.6, error=0.01, max_iter=100):
        """Usa Gauss-Siedel para resolver la red"""
        Vs = [[gen.V for gen in self.generadores if gen.bp == barra] for barra in self.barras]
        Vn = [vs[0] if vs else None for vs in Vs]
        Pn = [sum(carga.P for carga in self.cargas if carga.bp == barra) for barra in self.barras]
        Qn = [sum(carga.Q for carga in self.cargas if carga.bp == barra) for barra in self.barras]
        self._GS = GaussSiedel(self.barras, Vn, Pn, Qn, self.matrisAdmitancias, alpha)
        self._GS.resolver(error, max_iter)
        self._voltajes = self._GS.V[-1]
        self._potencias = [complex(p, q) if p and q else None for p, q in zip(self._GS.P, self._GS.Q[-1])]
        return self._GS


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


class AnalisisPerdidas:
    """Red simplificada de barras de despacho

    Determina las matrices de transformacion C
    """

    MAX_ERROR = 10 ** -4

    def __init__(self, red: RedPotencia):
        self.red = red
        self.barras = None
        self._C12 = None
        self._C23 = None
        self._C34 = None
        self._C4G = None
        self._C = None
        self._A = None
        self._B = None
        self._Ym = self.reduccionDeBarras()
        self.delta_x = np.zeros((1, len(self.red.generadores) + 1), dtype=np.complex64)

    @property
    def Ym(self):
        return self._Ym

    @property
    def Zm(self):
        self._zm = np.linalg.inv(self.Ym)
        return self._zm

    @property
    def Rb(self):
        self._Rb = self.Zm.real
        return self._Rb

    @property
    def C12(self):
        if self._C12 is None:
            self._C12 = self.separacionDeBarras()
        return self._C12

    @property
    def C23(self):
        if self._C23 is None:
            self._C23 = self.unificacionDeCargas()
        return self._C23

    @property
    def C34(self):
        if self._C34 is None:
            self._C34 = self.cambioDeReferencia()
        return self._C34

    @property
    def C4G(self):
        if self._C4G is None:
            self._C4G = self.cambioDeVariable()
        return self._C4G

    @property
    def C(self):
        if self._C is None:
            self._C = np.dot(np.dot(np.dot(self.C12, self.C23), self.C34), self.C4G)
        return self._C

    @property
    def A(self):
        self._A = np.dot(np.dot(self.C.transpose(), self.Rb), self.C)
        return self._A

    @property
    def B(self):
        return self.A.real

    def reduccionDeBarras(self):
        Ym = self.red.matrisAdmitancias
        self.barras = self.red.barras.copy()
        for barra in self.red.barrasTransmision:
            i = self.red.barras.index(barra)
            Ym = reduccionKron(Ym, i)
            self.barras = self.barras[:i] + self.barras[i + 1 :]
        return Ym

    def separacionDeBarras(self):
        """Separa las barras en barras de cargas y de generación"""
        cols = self.red.generadores + self.red.cargas
        C12 = np.zeros((len(self.barras), len(cols)), dtype=np.complex64)
        for y, barra in enumerate(self.barras):
            for x, elemento in enumerate(cols):
                if elemento.bp is barra:
                    C12[y][x] = 1
        return C12

    def unificacionDeCargas(self):
        """Unifica las cargas segun las corrientes que consumen"""
        ng = len(self.red.generadores)
        nc = len(self.red.cargas)
        if len(self.red.Is) != len(self.red.cargas):
            raise ValueError("Se requiere de una corriente por carga")
        k = [i / self.red.It for i in self.red.Is]
        C23 = np.zeros((ng + nc, ng + 1), dtype=np.complex64)
        for y in range(ng):
            C23[y][y] = 1
        for y, i in enumerate(k):
            C23[ng + y][ng] = i
        return C23

    def cambioDeReferencia(self):
        """Cambio de referencia de carga a tierra"""
        ng = len(self.red.generadores)
        C34 = np.identity(ng + 1, dtype=np.complex64)
        for x in range(ng):
            C34[ng][x] = -1
        return C34

    def cambioDeVariable(self):
        """Cambio de variable de corriente a potencia

        Args:
            IT (float): Corriente de tierra
            ms (list[float]): lista de factores de  trasnformacion
        """
        ng = len(self.red.generadores)
        C4G = np.identity(ng + 1, dtype=np.complex64)
        for n, (vg, ig) in enumerate(zip(self.red.Vg, self.red.Ig)):
            C4G[n][n] = (vg * ig).conjugate().real / ig
        C4G[ng][ng] = self.red.It
        return C4G

    def requiereIteracion(self):
        return max(abs(self.delta_x)) > self.MAX_ERROR


class Etapa:
    """Una etapa , combina un escenario con un grupo de combinaciones"""

    def __init__(self, nombre: str, duracion: float, escenario: Escenario, combinaciones: GrupoCombinaciones):
        self.nombre = nombre
        self.duracion = duracion
        self.escenario = escenario
        self.combinaciones = combinaciones
        self._costos = None

    def procesar(self, rp: RedPotencia):
        """Procesa la red de potencia en cada escenario y combinación, y obtiene su costo"""
        self._costos = {}
        for comb in self.combinaciones:
            print(f"Procesando Etapa {self}({self.escenario}, {comb}):")
            print()
            rp.procesar(comb, self.escenario)
            rp.resolver()
            print(rp.GaussSiedel)
            ap = AnalisisPerdidas(rp)
            print("Matris B")
            print(ap.B)
            print()
            de = DespachoEconomico(f"{self}-{self.escenario}-{comb}", comb)
            resultado = de.despacho(self.escenario)
            Ft = sum(gen["f"] for gen in resultado["despacho"])
            print(f"Resultado: Ft/$={Ft}")
            for gen in resultado["despacho"]:
                print("  {nombre}: P={p} , F/$={f}".format(**gen))
            self._costos[comb] = Ft * self.duracion
            print()

    @property
    def costos(self):
        if self._costos is None:
            self._costos = {comb: None for comb in self.combinaciones}
        return self._costos

    def __hash__(self):
        return self.nombre.__hash__()

    def __str__(self):
        return self.nombre

    def __repr__(self):
        return self.nombre


class GrupoEtapas(list):
    """Grupo de etapas"""

    def __init__(self, *etapas: Etapa):
        super().__init__(etapas)
        self.combinaciones = GrupoCombinaciones(*(comb for etapa in self for comb in etapa.combinaciones))

    def __str__(self):
        table = Texttable()
        table.set_deco(Texttable.HEADER | Texttable.VLINES)
        table.set_cols_align(["r"] + ["c"] * len(self))
        table.header(["Etapa"] + [etapa.nombre for etapa in self])
        table.add_row(["Escenario", *[etapa.escenario.nombre for etapa in self]])
        table.add_row(["Duración", *[etapa.duracion for etapa in self]])
        table.add_rows(
            [
                [comb.nombre, *[" ○ " if comb in etapa.combinaciones else " ❌ " for etapa in self]]
                for comb in self.combinaciones
            ],
            header=False,
        )
        return table.draw()

    def procesarCostos(self, rp: RedPotencia):
        for etapa in self:
            etapa.procesar(rp)
        table = Texttable()
        table.set_deco(Texttable.HEADER)
        table.set_cols_align(["r"] + ["c"] * len(self))
        table.header(["Etapa"] + [etapa.nombre for etapa in self])
        table.add_rows(
            [[comb.nombre, *[etapa.costos.get(comb, "❌") or "❌" for etapa in self]] for comb in self.combinaciones],
            header=False,
        )
        return table.draw()

    def determinarCaminoOptimo(self):
        """Determina el camino optimo usando el algoritmo de Dijkstra"""
        costos = {}
        previo = {}
        for etapa in self:
            for comb in etapa.combinaciones:
                nodo = (etapa, comb)
                costos[nodo] = inf
                previo[nodo] = None
        # analisis nodos primera etapa
        etapa = self[0]
        for comb in etapa.combinaciones:
            nodo = (etapa, comb)
            costos[nodo] = etapa.costos[comb]
        # analsis otros nodos
        for x, etapa_inicio in enumerate(self[:-1]):
            a = x + 1
            etapa_dest = self[a]
            for comb_inicio in etapa_inicio.combinaciones:
                nodo_inicio = (etapa_inicio, comb_inicio)
                for comb_dest in etapa_dest.combinaciones:
                    nodo_destino = (etapa_dest, comb_dest)
                    if (
                        etapa_dest.costos[comb_dest] is None
                        or self.combinaciones.transiciones[comb_inicio][comb_dest] is None
                    ):
                        # la combinacion no es factible
                        continue
                    costo = (
                        costos[nodo_inicio]
                        + self.combinaciones.transiciones[comb_inicio][comb_dest]
                        + etapa_dest.costos[comb_dest]
                    )
                    if costo <= costos[nodo_destino]:
                        costos[nodo_destino] = costo
                        previo[nodo_destino] = nodo_inicio
        # analisis nodos ultima etapa
        etapa = self[-1]
        comb = min(etapa.combinaciones, key=lambda comb: costos[(etapa, comb)])
        # mejor ruta
        ruta = []
        nodo = (etapa, comb)
        costo = costos[nodo]
        while nodo is not None:
            ruta.append(nodo)
            nodo = previo[nodo]
        return ruta, costo
