#!/usr/bin/env python3
# despacho.py
# Herramientas de Wolfang para calculos rapidos en sistemas de potencia 2

# Pip
from texttable import Texttable

# Local
from sistemasDePotencia.comun import paralelo, inf

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
class ElementoDespacho:
    """Elemento base

    Args:
        nombre (str): nombre para mostar del elemento
    """

    def __init__(self, nombre, **kwargs):
        self.nombre = nombre

    def __str__(self):
        return f"{self.nombre}"


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

    def __init__(self, nombre, *generadores, carga=None,  costo=1, **kwargs):
        self.gens = generadores
        self.costo = costo
        self._at = None
        self._bt = None
        super().__init__(nombre, **kwargs)

    def __str__(self):
        return (
            f"Barra {self.nombre}:\tλt = {self.at:0.4f}*Pt + {self.bt:0.2f},"
            f"\t{self.Pmax} >= Pt >= {self.Pmin}\n" + "".join(f"\t{g}\n" for g in self.gens)
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
            self._at = calc_a(*self.gens)
        return self._at

    @property
    def bt(self):
        if self._bt is None:
            self._bt = calc_b(self.at, *self.gens)
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


class Combinacion:
    def __init__(self, nombre, *generadores):
        self.nombre = nombre
        self.gens = generadores

    def __str__(self):
        return "Combinación {self.nombre}: " + " ".join(self.gens)


class GrupoCombinaciones:
    def __init__(self, *combinaciones):
        self.combs = combinaciones
        self._gens = None
        self._transiciones = None

    @property
    def gens(self):
        if self._gens is None:
            self._gens = list({gen for comb in self.combs for gen in comb.gens})
            self._gens.sort(key=lambda x: x.nombre)
        return self._gens

    def __str__(self):
        table = Texttable()
        table.set_deco(Texttable.HEADER)
        table.set_cols_align(["r"] + ["c"] * len(self.gens))
        table.add_rows(
            [
                ["Comb.", *(gen.nombre for gen in self.gens)],
                *([comb.nombre, *["✅" if gen in comb.gens else "❌" for gen in self.gens]] for comb in self.combs),
            ]
        )
        return table.draw()

    def _calc_costo_transicion(self, c1, c2):
        gens = {*c1.gens, *c2.gens}
        costo = 0
        for gen in gens:
            if gen not in c1.gens:
                costo += gen.Ccon
            elif gen not in c2.gens:
                costo += gen.Cdes
        return costo

    @property
    def transiciones(self):
        if self._transiciones is None:
            self._transiciones = []
            for c1 in self.combs:
                for c2 in self.combs:
                    if c1 is c2:
                        continue
                    self._transiciones.append((c1, c2, self._calc_costo_transicion(c1, c2)))
        return self._transiciones

    def reporteTransiciones(self):
        table = Texttable()
        table.set_deco(Texttable.HEADER | Texttable.VLINES)
        table.header(["Tran.", "Costo"])
        table.add_rows([[f"{c1.nombre} -> {c2.nombre}", costo] for c1, c2, costo in self.transiciones], header=False)
        return table.draw()



if __name__ == "__main__":
    test()
