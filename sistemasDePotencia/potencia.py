#!/usr/bin/env python3
# potencia.py
# Herramientas de Wolfang para calculos rapidos en sistemas de potencia 1

# Standard Library
from math import atan, sqrt, acos

# Pip
import numpy as np

from sistemasDePotencia.comun import (
    rect,
    degrees,
    cambio_base,
    display_rect,
    display_polar,
    display_single,
    S,
    V,
    I,
    Z,
    Y,
)


# Clases de ayuda


class DeltaEstrella:
    """Calcular rapidamente un Delta/Estrella en Impedancia o Admitancia

    Zabc: impedancias delta
    Yabc: admitancias delta
    Z123: impedancias estrella
    Y123: admitancias estrella
    1 2 3: nodos primario, secundario, terciario
    a b c: lineas entre nodos

    delta:

    Np---c---Ns
      \     /
       b   a
        \ /
         Nt

    strella:

    Np       Ns
       1   2
         .
         3
         Nt
    """

    Np = "Np"
    Ns = "Ns"
    Nt = "Nt"
    Za = None
    Zb = None
    Zc = None
    Z1 = None
    Z2 = None
    Z3 = None
    Ya = None
    Yb = None
    Yc = None
    Y1 = None
    Y2 = None
    Y3 = None

    def __init__(self, Np=None, Ns=None, Nt=None):
        if Np:
            self.Np = Np
        if Ns:
            self.Ns = Ns
        if Nt:
            self.Nt = Nt

    def _calc(self):
        """Usando las impedancia de estrella, calcula los otros valores"""
        self.Y1 = 1 / self.Z1
        self.Y2 = 1 / self.Z2
        self.Y3 = 1 / self.Z3
        _edi = (self.Z1 * self.Z2) + (self.Z2 * self.Z3) + (self.Z3 * self.Z1)
        self.Za = _edi / self.Z1
        self.Zb = _edi / self.Z2
        self.Zc = _edi / self.Z3
        _eda = self.Y1 + self.Y2 + self.Y3
        self.Ya = self.Y2 * self.Y3 / _eda
        self.Yb = self.Y3 * self.Y1 / _eda
        self.Yc = self.Y1 * self.Y2 / _eda

    def impedanciasEstrella(self, Z1, Z2, Z3):
        self.Z1 = Z1
        self.Z2 = Z2
        self.Z3 = Z3
        self._calc()

    def impedanciaDelta(self, Za, Zb, Zc):
        _dei = Za + Zb + Zc
        self.Z1 = Zb * Zc / _dei
        self.Z2 = Zc * Za / _dei
        self.Z3 = Za * Zb / _dei
        self._calc()

    def AdmitanciaEstrella(self, Y1, Y2, Y3):
        self.Z1 = 1 / Y1
        self.Z2 = 1 / Y2
        self.Z3 = 1 / Y3
        self._calc()

    def AdmitanciaDelta(self, Ya, Yb, Yc):
        _dea = (Ya * Yb) + (Yb * Yc) + (Yc * Ya)
        self.Z1 = Ya / _dea
        self.Z2 = Yb / _dea
        self.Z3 = Yc / _dea
        self._calc()

    def __str__(self):
        return (
            f"{self.Np}⟍      ⟋ {self.Ns}\n"
            f"  Z1⟍  ⟋ Z2    Z1: {display_rect(self.Z1)} Ω\n"
            f"      Y        Z2: {display_rect(self.Z2)} Ω\n"
            f"      |Z3      Z3: {display_rect(self.Z3)} Ω\n"
            f"      {self.Nt}\n"
            f"\n"
            f"{self.Np}---Zc---{self.Ns}\n"
            f"  \      /     Za: {display_rect(self.Za)} Ω\n"
            f"   Zb   Za     Zb: {display_rect(self.Zb)} Ω\n"
            f"    \  /       Zc: {display_rect(self.Zc)} Ω\n"
            f"     {self.Nt}\n"
            f"\n"
            f"{self.Np}⟍      ⟋ {self.Ns}\n"
            f"  Y1⟍  ⟋ Y2    Y1: {display_rect(self.Y1)} ℧\n"
            f"      Y        Y2: {display_rect(self.Y2)} ℧\n"
            f"      |Y3      Y3: {display_rect(self.Y3)} ℧\n"
            f"      {self.Nt}\n"
            f"\n"
            f"{self.Np}---Yc---{self.Ns}\n"
            f"  \      /     Ya: {display_rect(self.Ya)} ℧\n"
            f"   Yb   Ya     Yb: {display_rect(self.Yb)} ℧\n"
            f"    \  /       Yc: {display_rect(self.Yc)} ℧\n"
            f"     {self.Nt}\n"
        )


def ModeloPI(bp, bs, Yp, Yps, Ys):
    """Calcula el flujo de corriente y potencia en una serie de admitancias en forma de pi"""
    Yp0 = Y("Yp0", Yp, bp)
    Yps = Y("Yps", Yps, bp)
    Ys0 = Y("Ys0", Ys, bp)
    Ip0 = I("Ip0", bp.V * Yp0, bp)
    Is0 = I("Is0", bs.V * Ys0, bp)
    Ips = I("Ips", (bp.V - bs.V) * Yps, bp)
    Sps = S("Sps", bp.V * (Ip0 + Ips).conjugate(), bp)
    Ssp = S("Ssp", bs.V * (Is0 - Ips).conjugate(), bp)
    Sper = S("Sper", Sps + Ssp, bp)
    return (
        f"  {Yp0}\n"
        f"  {Yps}\n"
        f"  {Ys0}\n"
        f"  {Ip0}\n"
        f"  {Is0}\n"
        f"  {Ips}\n"
        f"  {Sps}\n"
        f"  {Ssp}\n"
        f"  {Sper}"
        f"""
    {bp.nombre}             {bs.nombre}
    P     → I →     S
      --+-[Yps]-+--
    ↓   |       |   ↓
    I [Yp0]   [Ys0] I
    ↓   |       |   ↓
      --+-------+--"""
    )


def tap_fic(Z, t):
    Y = 1 / Z
    return (
        f"  Z / t:       {display_rect(Z/t)} Ω pu\n"
        f"  Z / (1-t):   {display_rect(Z/(1-t))} Ω pu\n"
        f"  Z / t(t-1):  {display_rect(Z/(t**2-t))} Ω pu\n"
        f"  Y * t:       {display_rect(Y*t)} ℧ pu\n"
        f"  Y * (1-t):   {display_rect(Y*(1-t))} ℧ pu\n"
        f"  Z * t(t-1):  {display_rect(Y*(t**2-t))} ℧ pu\n"
    )


# Clases de elementos de potencia


class Elemento:
    """Elemento trifasico base

    Args:
        nombre: nombre unico para mostar del elemento
    """

    def __init__(self, nombre, **kwargs):
        self.nombre = nombre

    def __str__(self):
        return f"{self.nombre}"

    def __repr__(self):
        return f"{self.nombre}"

    def __hash__(self):
        return self.nombre.__hash__()


class Barra(Elemento):
    """Una Barra trifasica

    Establecen los valores bases para los elementos conectados a ellos

    Args:
        Sb: Potencia base
        Vb: Voltaje Base
    """

    def __init__(self, nombre, *, Sb, Vb, Ib=None, Zb=None, Yb=None, V=None):
        self.Sb = Sb
        self.Vb = Vb
        self.Ib = Ib or self.Sb / self.Vb / sqrt(3)
        self.Zb = Zb or self.Vb ** 2 / self.Sb
        self.Yb = Yb or 1 / self.Zb
        self.V = V
        super().__init__(nombre)

    def __str__(self):
        return (
            f"{self.nombre}:\n"
            f"  Sb: {display_single(self.Sb)} W\n"
            f"  Vb: {display_single(self.Vb)} V\n"
            f"  Ib: {display_single(self.Ib)} A\n"
            f"  Zb: {display_single(self.Zb)} Ω\n"
            f"  Yb: {display_single(self.Yb)} ℧\n"
        )


class GeneradorIdeal(Elemento):
    """Generador Ideal conectado de tierra a barra primaria"""

    def __init__(self, nombre, bp, *, Pn, Vn, Qnmin=None, Qnmax=None, **kwargs):
        super().__init__(nombre, **kwargs)
        self.bp = bp
        self.Pn = S("Pn", Pn / bp.Sb, bp)
        self.Vn = V("Vn", Vn / bp.Vb, bp)
        self.Qnmin = None if Qnmin is None else S("Qmin", complex(0, Qnmin) / bp.Sb, bp)
        self.Qnmax = None if Qnmax is None else S("Qmax", complex(0, Qnmax) / bp.Sb, bp)

    def __str__(self):
        return (
            f"{self.nombre}: ({self.bp.nombre})\n"
            f"  {self.Pn}\n"
            f"  {self.Vn}\n"
            f" Valores:\n"
            f"  {self.Qnmin}\n"
            f"  {self.Qnmax}\n"
        )


class GeneradorBasico(GeneradorIdeal):
    """Generador conectado a la barra bp con parametro en absoluto"""

    def __init__(self, nombre, *, bp, Sn, PF, Vn, R=None, X=None, Zev=None):
        self.PF = PF
        Sc = rect(Sn, acos(PF))
        self.S = S("Sn", Sc / bp.Sb, bp)
        super().__init__(
            nombre,
            bp,
            Pn=Sc.real,
            Vn=Vn,
        )
        self.bp = bp
        self.Vn = Vn
        self.V = V("V", self.Vn / self.bp.Vb, self.bp)
        self.R = Z("R", complex(R, 0), self.bp)
        self.X = Z("X", complex(0, X), self.bp)
        self.Z = Z("Z", self.R + self.X, self.bp)
        self.Y = Y("Y", 1 / self.Z, self.bp)

    def __str__(self):
        return (
            f"{self.nombre}: ({self.bp.nombre})\n"
            f"  {self.S}\n"
            f"  {self.V}\n"
            f" Valores:\n"
            f"  {self.Z}\n"
            f" Modelo Admitancias:\n"
            f"  {self.Y}\n"
        )


class Linea(Elemento):
    """Una Linea Mediana de transmicion trifasica

    Interconecta dos barras con los mismos valores bases

    Args:
        bp: Barra del lado primario
        bp: Barra del lado secundario
        R: Paramtro resistivo en serie
        X: Paramtro inductivo en serie
        G: Paramtro conductivo en paralelo
        B: Paramtro suceptivo en paralelo
    """

    def __init__(self, nombre, *, bp, bs, R, X, G, B):
        super().__init__(nombre)
        self.bp = bp
        self.bs = bs
        self.R = Z("R", complex(R, 0), self.bp)
        self.X = Z("X", complex(0, X), self.bp)
        self.G = Y("G", complex(G, 0), self.bp)
        self.B = Y("B", complex(0, B), self.bp)
        self.Yps = Y("Yps", 1 / (self.R + self.X), self.bp)
        self.Yp0 = Y("Yp0", (self.G + self.B) / 2, self.bp)
        self.Ys0 = Y("Ys0", (self.G + self.B) / 2, self.bp)

    def __str__(self):
        return (
            f"{self.nombre}: ({self.bp.nombre} → {self.bs.nombre})\n"
            f" Valores:\n"
            f"  {self.R}\n"
            f"  {self.X}\n"
            f"  {self.G}\n"
            f"  {self.B}\n"
            f" Modelo Admitancias:\n"
            f"  {self.Yps}\n"
            f"  {self.Yp0}\n"
            f"  {self.Ys0}\n"
        )

    def en_real(self):
        return (
            f"{self.nombre}:\n"
            f" Valores Reales:\n"
            f"  {self.R.en_real()}\n"
            f"  {self.X.en_real()}\n"
            f"  {self.G.en_real()}\n"
            f"  {self.B.en_real()}\n"
            f" Modelo Admitancias en Real:\n"
            f"  {self.bp.nombre} -> {self.bs.nombre}: {self.Yps.en_real()}\n"
            f"  {self.bp.nombre} -> 0: {self.Yp0.en_real()}\n"
            f"  {self.bs.nombre} -> 0: {self.Ys0.en_real()}\n"
        )

    def modelo_admitancias(self):
        return f"Flujo de Potencias en {self.nombre}:\n" f"{ModeloPI(self.bp, self.bs, self.Yp0, self.Yps, self.Ys0)}\n"

    def flujo_potencias(self):
        return f"Flujo de Potencias en {self.nombre}:\n" f"{ModeloPI(self.bp, self.bs, self.Yp0, self.Yps, self.Ys0)}\n"


class LineaReal(Linea):
    """Una Linea Mediana de transmicion trifasica

    Usa los valores absolutos por kilometro

    Args:
        Sn: Potencia Nominal
        Vn: Voltaje Nominal
        In: Corriente Nominal
        L: Longitud (Kilometros)
        T: Numero de cables
        Rn: Parametro de Resistencia de la linea en valor absoluto (por kilmetro)
        Xn: Parametro de Reactancia de la linea en valor absoluto (por kilmetro)
        Gn: Parametro de Conductancia de la linea en valor absoluto (por kilmetro)
        Bn: Parametro de Suceptancia de la linea en valor absoluto (por kilmetro)
    """

    def __init__(self, nombre, *, Rn, Xn, Gn, Bn, bp, bs, T=1, L=1, Sn=None, Vn=None, In=None):
        self.Sn = Sn
        self.Vn = Vn
        self.In = In
        self.L = L
        self.T = T
        self.Rn = Rn
        self.Xn = Xn
        self.Gn = Gn
        self.Bn = Bn

        super().__init__(
            nombre,
            bp=bp,
            bs=bs,
            R=Rn * L / T / bp.Zb,
            X=Xn * L / T / bp.Zb,
            G=Gn * L * T / bp.Yb,
            B=Bn * L * T / bp.Yb,
        )

    def __str__(self):
        return (
            f"{self.nombre}: ({self.bp.nombre} → {self.bs.nombre}, {self.L} km, {self.T} Ternas)\n"
            f" Valores:\n"
            f"  {self.R}\n"
            f"  {self.X}\n"
            f"  {self.G}\n"
            f"  {self.B}\n"
            f" Modelo Admitancias:\n"
            f"  {self.Yps}\n"
            f"  {self.Yp0}\n"
            f"  {self.Ys0}\n"
        )


class BancoCondensadores(Elemento):
    def __init__(self, nombre, *, Qn, In, Prel, bp, bs):
        super().__init__(nombre)
        self.Qn = Qn
        self.In = In
        self.Prel = Prel
        self.bp = bp
        self.bs = bs
        self.Pn = self.Prel * self.Qn / sqrt(1 - self.Prel ** 2)
        self.Sn = complex(self.Pn, self.Qn)
        self.Vn = self.Sn / sqrt(3) / self.In
        self.Zn = self.Vn ** 2 / self.Sn
        self.Z = self.Zn / self.bp.Zb
        self.Y = 1 / self.Z

    def __str__(self):
        return (
            f"{self.nombre}:\n"
            f"  Sn: {display_rect(self.Sn)} VA\n"
            f"  Zn: {display_rect(self.Zn)} Ω\n"
            f"  Vn: {display_polar(self.Vn)} V\n"
            f"  Z:  {display_rect(self.Z)} Ω pu\n"
            f"  Y:  {display_rect(self.Y)} ℧ pu\n"
        )


class TransformadorSimple(Elemento):
    """Transformador Trifasico simple de 2 devanados

    tiene un lado primario y un lado secundario,
    no hay preferencia de que uno sea Alta o Baja,
    eso depende de las barra a la que estan conectados

    Args:
        Sn: Potencia Nominal
        Vp: Voltaje nominal del lado primario
        Vs: Voltaje nominal del lado secundario

        Zev: Impedancia equivalente en pu
        rel: relacion X/R

        bp: Barra del lado primario
        bs: Barra del lado secundario

    """

    def __init__(self, nombre, *, Sn, Vp, Vs, Zev, rel, bp, bs):
        super().__init__(nombre)
        self.Sn = Sn
        self.Vp = Vp
        self.Vs = Vs
        self.Zev = Z("Zev", rect(Zev, degrees(atan(rel))))
        self.Zev.muestra = "polar"
        self.bp = bp
        self.bs = bs
        self.Zen = Z(
            "Zen",
            self.Zev * ((self.Vp / self.bp.Vb) ** 2) * (self.bp.Sb / self.Sn),
            self.bp,
        )
        self.Yen = Y("Yen", 1 / self.Zen, self.bp)

    def __str__(self):
        return (
            f"{self.nombre}: ({self.bp.nombre} → {self.bs.nombre})\n"
            f"  Sn:  {display_single(self.Sn)} VA\n"
            f"  Vn:    {display_single(self.Vp)} / {display_single(self.Vs)} V\n"
            f"  Bases: {self.bp.nombre:>11s} / {self.bs.nombre:>11s}\n"
            f" Valores:\n"
            f"  {self.Zev}\n"
            f"  {self.Zen}\n"
            f"  {self.Yen}\n"
        )


class TransformadorTap(TransformadorSimple):
    """Transformador con tap, basado en Transformador Simple

    Se asume que el tap esta conectado al primario

    Args:
        dt: regulacion por cada paso, Ej: 20% y 32 pasos -> dt = 0.2 / 32 = 0.00625
        pos: paso actual, 0 es el tap nominal, +1 es el tap por ensima, -1 el tap por devajo
    """

    def __init__(self, nombre, *, dt, pos, **kwargs):
        self.pos = pos
        self.dt = dt
        super().__init__(nombre, **kwargs)
        del self.Zen, self.Yen
        self.t = 1 + (self.dt * self.pos)
        self.Zp = Z(
            "Zp",
            self.Zev * ((self.Vp * self.t / self.bp.Vb) ** 2) * (self.bp.Sb / self.Sn),
            self.bp,
        )
        self.Zs = Z(
            "Zs",
            self.Zev * ((self.Vs / self.bs.Vb) ** 2) * (self.bs.Sb / self.Sn),
            self.bs,
        )
        self.Yp = Y("Yp", 1 / self.Zp, self.bp)
        self.Ys = Y("Ys", 1 / self.Zs, self.bs)
        self.Yps = Y("Yps", self.Yp * self.t, self.bp)
        self.Ys0 = Y("Ys0", self.Yp * (1 - self.t), self.bp)
        self.Yp0 = Y("Yp0", self.Yp * (self.t ** 2 - self.t), self.bp)

    def __str__(self):
        return (
            f"{self.nombre}: ({self.bp.nombre} → {self.bs.nombre})\n"
            f"  Sn: {display_single(self.Sn)} VA\n"
            f"  Vn: {display_single(self.Vp)} / {display_single(self.Vs)} V\n"
            f"  t:  {self.t} (Tap {'+' if self.pos >=0 else '-'}{abs(self.pos)})\n"
            f" Valores:\n"
            f"  {self.Zev}\n"
            f"  {self.Zp}\n"
            f"  {self.Zs}\n"
            f"  {self.Yp}\n"
            f"  {self.Ys}\n"
            f" Modelo Admitancias:\n"
            f"  {self.Yps}\n"
            f"  {self.Yp0}\n"
            f"  {self.Ys0}\n"
        )

    def flujo_potencias(self):
        return f"Flujo de Potencias en {self.nombre}:\n" f"{ModeloPI(self.bp, self.bs, self.Yp0, self.Yps, self.Ys0)}\n"


class TransformadorTriple(Elemento):
    """Transformador Trifasico de 3 devanados

    Tiene un lado primario, secundario y terciario

    Args:

    """

    def __init__(
        self,
        nombre,
        *,
        Snp,
        Sns,
        Snt,
        Vnp,
        Vns,
        Vnt,
        bp,
        bs,
        bt,
        Zept=None,
        Zeps=None,
        Zest=None,
    ):
        super().__init__(nombre)
        self.Snp = Snp
        self.Sns = Sns
        self.Snt = Snt
        self.Vnp = Vnp
        self.Vns = Vns
        self.Vnt = Vnt
        self.bp = bp
        self.bs = bs
        self.bt = bt
        self.Zept = Zept
        self.Zeps = Zeps
        self.Zest = Zest
        self.Zp = None
        self.Zs = None
        self.Zt = None
        self.Yp = None
        self.Ys = None
        self.Yt = None
        self.Zpt = self.Zept and cambio_base(self.Zept, self.Snp, self.Vnp, self.bp.Sb, self.bp.Vb)
        self.Zps = self.Zeps and cambio_base(self.Zeps, self.Snp, self.Vns, self.bs.Sb, self.bs.Vb)
        self.Zst = self.Zest and cambio_base(self.Zest, self.Snp, self.Vnt, self.bt.Sb, self.bt.Vb)
        self.Ypt = self.Zpt and 1 / self.Zpt
        self.Yps = self.Zps and 1 / self.Zps
        self.Yst = self.Zst and 1 / self.Zst
        self.tps = (self.bp.Vb / self.Vnp) / (self.bs.Vb / self.Vns)
        self.tpt = (self.bp.Vb / self.Vnp) / (self.bt.Vb / self.Vnt)
        self.tst = (self.bs.Vb / self.Vns) / (self.bt.Vb / self.Vnt)

    def pruebas_corto_circuito(self, *, Vpt, Ipt, Vps, Ips, Vst, Ist):
        zpt = (Vpt / self.bp.Vb) / (Ipt / self.bt.Ib)
        zps = (Vps / self.bp.Vb) / (Ips / self.bs.Ib)
        zst = (Vst / self.bs.Vb) / (Ist / self.bt.Ib)
        sumas = np.array([[zpt], [zps], [zst]])
        matris = np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1]])
        resultados = np.linalg.inv(matris).dot(sumas)
        [Zp], [Zs], [Zt] = resultados
        self.Zp = rect(Zp, 90)
        self.Zs = rect(Zs, 90)
        self.Zt = rect(Zt, 90)
        self.Yp = 1 / self.Zp
        self.Ys = 1 / self.Zs
        self.Yt = 1 / self.Zt

    def __str__(self):
        return (
            (
                f"{self.nombre}:\n"
                f"  Sn:    {display_single(self.Snp)} / {display_single(self.Sns)} / {display_single(self.Snt)} VA\n"
                f"  Vn:    {display_single(self.Vnp)} / {display_single(self.Vns)} / {display_single(self.Vnt)} V\n"
                f"  Bases: {self.bp.nombre:>11s} / {self.bs.nombre:>11s} / {self.bt.nombre:>11s}\n"
                f"  tps:   {self.tps}\n"
                f"  tpt:   {self.tpt}\n"
                f"  tst:   {self.tst}\n"
            )
            + (f"  Zpt: {display_rect(self.Zpt)} Ω pu\n" if self.Zpt else "")
            + (f"  Zps: {display_rect(self.Zps)} Ω pu\n" if self.Zps else "")
            + (f"  Zst: {display_rect(self.Zst)} Ω pu\n" if self.Zst else "")
            + (f"  Ypt: {display_rect(self.Ypt)} Ω pu\n" if self.Ypt else "")
            + (f"  Yps: {display_rect(self.Yps)} Ω pu\n" if self.Yps else "")
            + (f"  Yst: {display_rect(self.Yst)} Ω pu\n" if self.Yst else "")
            + (f"  Zp: {display_rect(self.Zp)} Ω pu\n" if self.Zp else "")
            + (f"  Zs: {display_rect(self.Zs)} Ω pu\n" if self.Zs else "")
            + (f"  Zt: {display_rect(self.Zt)} Ω pu\n" if self.Zt else "")
            + (f"  Yp: {display_rect(self.Yp)} Ω pu\n" if self.Yp else "")
            + (f"  Ys: {display_rect(self.Ys)} Ω pu\n" if self.Ys else "")
            + (f"  Yt: {display_rect(self.Yt)} Ω pu\n" if self.Yt else "")
            + (f"  Zp + Zs: {display_rect(self.Zp+self.Zs)} Ω pu\n" if self.Zp else "")
            + (f"  Zp + Zt: {display_rect(self.Zp+self.Zt)} Ω pu\n" if self.Zs else "")
            + (f"  Zs + Zt: {display_rect(self.Zt+self.Zs)} Ω pu\n" if self.Zt else "")
        )


class BancoCapasitores(Elemento):
    """Banco de Capasitores conectados a Tierra"""

    def __init__(self, nombre, *, bp, Sn, Vn, pj, **kwargs):
        self.bp = bp
        self.Sn = S("Sn", Sn + abs(Sn) * pj, bp)
        self.Vn = V("Vn", Vn, self.bp)
        self.In = I("In", self.Sn / sqrt(3) / self.Vn, self.bp)
        self.Zn = Z("Zn", self.Vn ** 2 / self.Sn, self.bp)
        self.Yn = Y("Yn", 1 / self.Zn, self.bp)
        super().__init__(nombre, **kwargs)

    def __str__(self):
        return (
            f"{self.nombre}: ({self.bp.nombre})\n"
            f"  {self.Sn}\n"
            f"  {self.Vn}\n"
            f"  {self.In}\n"
            f"  {self.Zn}\n"
            f"  {self.Yn}\n"
        )

    def en_real(self):
        return (
            f"{self.nombre}: ({self.bp.nombre})\n"
            f"  {self.Sn.en_real()}\n"
            f"  {self.Vn.en_real()}\n"
            f"  {self.In.en_real()}\n"
            f"  {self.Zn.en_real()}\n"
            f"  {self.Yn.en_real()}\n"
        )


class CargaIdeal(Elemento):
    """Carga Ideal conectado de tierra a barra primaria"""

    def __init__(self, nombre, *, bp, Pn, Qn, **kwargs):
        self.bp = bp
        self.P = S("P", Pn / bp.Sb, bp)
        self.Q = S("Q", Qn * 1j / bp.Sb, bp)
        self.S = S("S", complex(Pn, Qn) / bp.Sb, bp)
        self.Y = Y("Y", 0, bp)
        self.S.muestra = "polar"
        super().__init__(nombre, **kwargs)

    def __str__(self):
        return f"{self.nombre}: ({self.bp.nombre})\n" f" Valores:\n" f"  {self.S}\n" f"  {self.P}\n" f"  {self.Q}\n"
