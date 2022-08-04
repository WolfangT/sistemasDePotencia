#!/usr/bin/env python3
# despacho.py
# Herramientas de Wolfang para calculos rapidos en sistemas de potencia 2

from math import sin, acos

# Pip
import numpy as np
from texttable import Texttable

from sistemasDePotencia.comun import inf, paralelo
from sistemasDePotencia.potencia import Elemento

# Funciones de ayuda


def calc_a(*gens):
    """Calcula el parametro a total de una serie de generadores"""
    return paralelo(gen.a for gen in gens)


def calc_b(a, *gens):
    """Calcula el parametro b total de una serie de generadores"""
    return a * sum(gen.b / gen.a for gen in gens)


def tabla_partes_iguales(data):
    """Dibuja una tabla de Partes iguales a base de un analisis"""
    numero_gen = len(data["gens"])
    table = Texttable()
    table.set_deco(Texttable.HEADER)
    table.set_cols_dtype(["i"] + ["f"] * numero_gen + ["i"] * numero_gen)
    table.add_rows(
        [
            ["Pt"] + [f"λ{i.nombre}" for i in data["gens"]] + [f"P{i.nombre}" for i in data["gens"]],
            *(
                [pt, *(g["lambda"] for g in pc), *(g["p"] for g in pc)]
                for pt, pc in zip(data["potencias"], data["partes_iguales"])
            ),
        ]
    )
    return table.draw()


def tabla_despacho_economnico(data):
    """Dibuja una tabla de Despacho Economico a base de un analisis"""
    # numero_gen = len(data["gens"])
    table = Texttable()
    table.set_deco(Texttable.HEADER)
    table.add_rows(
        [
            ["Pt"] + [f"λ{i.nombre}" for i in data["gens"]] + [f"P{i.nombre}" for i in data["gens"]],
            *(
                [pt, *(g["lambda"] for g in pc), *(g["p"] for g in pc)]
                for pt, pc in zip(data["potencias"], data["despacho"])
            ),
        ]
    )
    return table.draw()


def tabla_costos(data):
    """Dibuja una tabla de Costos a base de un analisis"""
    numero_gen = len(data["gens"])
    table = Texttable()
    table.set_deco(Texttable.HEADER)
    table.set_cols_dtype(["i"] + ["f"] * numero_gen + ["f", "f"])
    table.add_rows(
        [
            ["Pt"] + [f"Δf{i.nombre}" for i in data["gens"]] + ["ΔF", "Ahorro Anual"],
            *(
                [pt, *(g for g in pc["delta_f"]), pc["delta_ft"], pc["ahorro"]]
                for pt, pc in zip(data["potencias"], data["costos"])
            ),
        ]
    )
    return table.draw()


# Clases

# Elementos de Red
class ElementoDespacho(Elemento):
    """Elemento base par adespacho"""


class CargaDespacho(ElementoDespacho):
    """Representa una carga, sus caracteristicas de consumo dependen del escenario"""

    def __str__(self):
        return f"Carga {self.nombre}"


class GeneradorDespacho(ElementoDespacho):
    """Representa a un generador

    Args:
        nombre (str): nombre para mostar del elemento
        a (int): parametro a del costo incremental
        b (int): parametro b del costo incremental
        c (int, optional): parametro c del costo incremental. Defaults to 0.
        Pmin (int): potencia minima a la que opera el generador. Defaults to 0.
        Pmax (int): potencia minima a la que opera el generador. Detaults to +inf
        Ccon (int): Costo de coneccion. Defaults to 0.
        Cdes (int): Costo de desconeccion. Defaults to 0.
    """

    def __init__(self, nombre, *, a, b, c=0, Pmin=0, Pmax=inf, Ccon=0, Cdes=0, **kwargs):
        self.a = a
        self.b = b
        self.c = c
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.Ccon = Ccon
        self.Cdes = Cdes
        super().__init__(nombre, **kwargs)

    def fijar_potencia(self, p):
        if p > self.Pmax:
            p = self.Pmax
        elif p < self.Pmin:
            p = self.Pmin
        lambda_ = self.a * p + self.b
        f = self.a / 2 * p ** 2 + self.b * p + self.c
        return {"lambda": lambda_, "p": p, "f": f, "nombre": self.nombre}

    def fijar_lambda(self, lambda_):
        p = (lambda_ - self.b) / self.a
        if p > self.Pmax:
            p = self.Pmax
            lambda_ = self.a * p + self.b
        elif p < self.Pmin:
            p = self.Pmin
            lambda_ = self.a * p + self.b
        f = self.a / 2 * p ** 2 + self.b * p + self.c
        return {"lambda": lambda_, "p": p, "f": f, "nombre": self.nombre}

    def fijar_0(self):
        return {"lambda": "-", "p": "-", "f": "-", "nombre": self.nombre}

    def __str__(self):
        return f"Gen {self.nombre}:\tλ = {self.a:0.4f}*P{self.nombre} + {self.b:5.2f},\t{self.Pmax} >= P{self.nombre} >= {self.Pmin}"


class BarraDespacho(ElementoDespacho):
    """Barra de generadores

    Args:
        generadores (list[Generador]): lista de generadores conectados
        carga (Carga): Carga conectada a la barra
        costo: Costo por unidad de combustible en los generadores de esta barra
    """

    def __init__(self, nombre, *generadores, carga=None, costo=1, **kwargs):
        self.gens = generadores
        self.costo = costo
        self.carga = carga
        self._at = None
        self._bt = None
        super().__init__(nombre, **kwargs)

    def __str__(self):
        return (
            f"Barra {self.nombre}:\tλt = {self.at:0.4f}*Pt + {self.bt:0.2f},"
            f"\t{self.Pmax} >= Pt >= {self.Pmin}\n"
            + "".join(f"\t{g}\n" for g in self.gens)
            + (f"\t{self.carga}\n" if self.carga else "")
        )

    def _ordenar(self, data):
        """Obtiene el index de un resultado segun la lista interna de generadores"""
        for gen in self.gens:
            if gen.nombre == data["nombre"]:
                return self.gens.index(gen)

    @property
    def Pmin(self):
        return sum(gen.Pmin for gen in self.gens)

    @property
    def Pmax(self):
        return sum(gen.Pmax for gen in self.gens)

    @property
    def at(self):
        if self._at is None:
            self._at = calc_a(*self.gens) if self.gens else inf
        return self._at

    @property
    def bt(self):
        if self._bt is None:
            self._bt = calc_b(self.at, *self.gens) if self.gens else inf
        return self._bt

    def calc_partes_iguales(self, pt, gens=None, fijo=None):
        """Calculo de generadores para la potencia pt en partes iguales"""
        if gens is None:
            gens = list(self.gens)
        if fijo is None:
            fijo = []
        p_pt = pt / len(gens)
        partes_iguales = []
        for gen in gens.copy():
            res = gen.fijar_potencia(p_pt)
            if res["p"] != p_pt:
                gens.remove(gen)
                pt -= res["p"]
                fijo.append(res)
                break
            partes_iguales.append(res)
        else:
            resultados = fijo + partes_iguales
            resultados.sort(key=self._ordenar)
            return resultados
        return self.calc_partes_iguales(pt, gens, fijo)

    def calc_despacho_economico(self, pt, gens=None, fijo=None):
        """Calculo de generadores para la potencia pt en despacho economico"""
        recalc = False
        if gens is None:
            gens = list(self.gens)
        if fijo is None:
            fijo = []
        a = calc_a(*gens)
        b = calc_b(a, *gens)
        lambda_t = a * pt + b
        despacho_economico = []
        for gen in gens.copy():
            res = gen.fijar_lambda(lambda_t)
            if res["lambda"] != lambda_t:
                recalc = True
            despacho_economico.append(res)
        if recalc:
            despacho_economico.sort(key=lambda x: x["lambda"], reverse=True)
            res = despacho_economico[0]
            fijo.append(res)
            for gen in gens:
                if gen.nombre == res["nombre"]:
                    break
            gens.remove(gen)
            pt -= res["p"]
            return self.calc_despacho_economico(pt, gens, fijo)
        else:
            resultados = fijo + despacho_economico
            resultados.sort(key=self._ordenar)
            return resultados

    def despacho(self, pt):
        """Crea el resultado de despachar la potencia pt tanto por partes
        iguales como por despacho economico
        """
        if pt < self.Pmin or pt > self.Pmax:
            return {
                "pt": pt,
                "partes_iguales": [gen.fijar_0() for gen in self.gens],
                "despacho": [gen.fijar_0() for gen in self.gens],
                "costos": {
                    "delta_f": ["-" for gen in self.gens],
                    "delta_ft": "-",
                    "ahorro": "-",
                },
            }
        # Partes iguales
        partes_iguales = self.calc_partes_iguales(pt)
        # Despacho
        despacho = self.calc_despacho_economico(pt)
        # Analisis de costos
        delta_f = []
        for gen_igual, gen_despacho in zip(partes_iguales, despacho):
            delta_f.append(gen_despacho["f"] - gen_igual["f"])
        delta_ft = sum(delta_f)
        ahorro = -self.costo * delta_ft
        analisis = {"delta_f": delta_f, "delta_ft": delta_ft, "ahorro": ahorro}
        return {
            "pt": pt,
            "partes_iguales": partes_iguales,
            "despacho": despacho,
            "costos": analisis,
        }

    def analisis(self, *p):
        """Crea la data analisando el despacho de la barra en diferentes potencias"""
        data = {
            "barra": self.nombre,
            "gens": self.gens,
            "potencias": [],
            "partes_iguales": [],
            "despacho": [],
            "costos": [],
        }
        for pt in p:
            calc = self.despacho(pt)
            data["potencias"].append(pt)
            data["partes_iguales"].append(calc["partes_iguales"])
            data["despacho"].append(calc["despacho"])
            data["costos"].append(calc["costos"])
        return data


# Combinaciones de elementos


class Combinacion(list):
    """Combinaciones de generadores de Despacho"""

    def __init__(self, nombre, *generadores):
        self.nombre = nombre
        super().__init__(generadores)

    @property
    def Pmax(self):
        return sum(gen.Pmax for gen in self)

    @property
    def Pmin(self):
        return sum(gen.Pmin for gen in self)

    def __str__(self):
        return self.nombre

    def __repr__(self):
        return self.nombre

    def __hash__(self):
        return self.nombre.__hash__()


class GrupoCombinaciones(list):
    """Grupo de combinaciones de generadores de despacho

    Calcula las transiciones entre ellas y sus costos
    """

    def __init__(self, *combinaciones):
        combs = list(set(combinaciones))
        combs.sort(key=lambda x: x.nombre)
        super().__init__(combs)
        self._gens = None
        self._transiciones = None

    @property
    def gens(self):
        """Lista de generadores en todas la combinaciones"""
        if self._gens is None:
            self._gens = list({gen for comb in self for gen in comb})
            self._gens.sort(key=lambda x: x.nombre)
        return self._gens

    def __str__(self):
        table = Texttable()
        table.set_deco(Texttable.HEADER)
        table.set_cols_align(["r"] + ["c"] * len(self.gens) + ["c", "c"])
        table.add_rows(
            [
                ["Comb.", *(gen.nombre for gen in self.gens), "Pmin", "Pmax"],
                *(
                    [comb.nombre, *["✅" if gen in comb else "❌" for gen in self.gens], comb.Pmin, comb.Pmax]
                    for comb in self
                ),
            ]
        )
        return table.draw()

    def _calc_costo_transicion(self, c1, c2):
        gens = {*c1, *c2}
        costo = 0
        for gen in gens:
            if gen not in c1:
                costo += gen.Ccon
            elif gen not in c2:
                costo += gen.Cdes
        return costo

    @property
    def transiciones(self):
        if self._transiciones is None:
            self._transiciones = {}
            for c1 in self:
                self._transiciones[c1] = {}
                for c2 in self:
                    self._transiciones[c1][c2] = self._calc_costo_transicion(c1, c2)
        return self._transiciones

    def reporteTransiciones(self):
        table = Texttable()
        table.set_deco(Texttable.HLINES | Texttable.VLINES | Texttable.HEADER)
        table.header(["Transición.", *[str(comb) for comb in self]])
        table.add_rows(
            [[str(c1), *[self.transiciones[c1][c2] for c2 in self.transiciones[c1]]] for c1 in self.transiciones],
            header=False,
        )
        return table.draw()


class Escenario(dict):
    """Escenario de cargas"""

    def __init__(self, nombre: str, cargas: dict[CargaDespacho, list[float, float]]):
        self.nombre = nombre
        super().__init__({carga: complex(s * pf, round(s * sin(acos(pf)), 2)) for carga, (s, pf) in cargas.items()})
        self.total = sum(s for s in self.values())

    def __str__(self):
        return (
            f"Escenario {self.nombre}:\n"
            + "\n".join(f"  {carga.nombre}: {s}" for carga, s in self.items())
            + "\n"
            + f" {self.total}"
        )

    def __hash__(self):
        return self.nombre.__hash__()


class GrupoEscenarios(list):
    """Grupo de Escenarios de Carga"""

    def __init__(self, *escenarios: Escenario):
        super().__init__(escenarios)
        self.cargas = {carga for es in self for carga in es}

    def __str__(self):
        table = Texttable()
        table.set_deco(Texttable.HEADER)
        table.set_cols_align(["r"] + ["c"] * len(self) + ["l"])
        table.header(["Escenario"] + [carga.nombre for carga in self.cargas] + ["Total"])
        table.add_rows(
            [[es.nombre, *[es.get(carga, "❌") for carga in self.cargas], es.total] for es in self],
            header=False,
        )
        return table.draw()


class Etapa:
    """Una etapa , combina un escenario con un grupo de combinaciones"""

    def __init__(self, nombre: str, duracion: float, escenario: Escenario, combinaciones: GrupoCombinaciones):
        self.nombre = nombre
        self.duracion = duracion
        self.escenario = escenario
        self.combinaciones = combinaciones
        self._costos = None

    def procesar(self):
        pass

    @property
    def costos(self):
        if self._costos is None:
            self._costos = {comb: 1000 * self.duracion for comb in self.combinaciones}
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
        table.set_deco(Texttable.HEADER)
        table.set_cols_align(["r"] + ["c"] * len(self))
        table.header(["Etapa"] + [etapa.nombre for etapa in self])
        table.add_row(["Escenario", *[etapa.escenario.nombre for etapa in self]])
        table.add_row(["Duración", *[etapa.duracion for etapa in self]])
        table.add_rows(
            [
                [comb.nombre, *["☆" if comb in etapa.combinaciones else "❌" for etapa in self]]
                for comb in self.combinaciones
            ],
            header=False,
        )
        return table.draw()

    def procesarCostos(self):
        for etapa in self:
            etapa.procesar()
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
