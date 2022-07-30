# Pip
import numpy as np

from sistemasDePotencia.comun import M, Vb13, Vb24, u, inf
from sistemasDePotencia.redes import (
    Barra,
    Linea,
    LineaReal,
    CargaIdeal,
    GaussSiedel,
    GeneradorIdeal,
    GeneradorBasico,
    BancoCapasitores,
    TransformadorTap,
    reduccionKron,
)
from sistemasDePotencia.despacho import (
    Combinacion,
    BarraDespacho,
    CargaDespacho,
    GeneradorDespacho,
    GrupoCombinaciones,
    tabla_costos,
    tabla_partes_iguales,
    tabla_despacho_economnico,
)


def trabajo_final():
    Sb = 1000 * M
    SbTxA = 1500 * M
    SbTxB = 500 * M
    print(B1 := Barra("B1", Sb=Sb, Vb=Vb13))
    print(B2 := Barra("B2", Sb=Sb, Vb=Vb24))
    print(B3 := Barra("B3", Sb=Sb, Vb=Vb24))
    print(B4 := Barra("B4", Sb=Sb, Vb=Vb13))
    print(B5 := Barra("B5", Sb=Sb, Vb=Vb24))
    print(B6 := Barra("B6", Sb=Sb, Vb=Vb13))
    print(B7 := Barra("B7", Sb=Sb, Vb=Vb24))
    print(B8 := Barra("B8", Sb=Sb, Vb=Vb13))
    print(G1 := GeneradorBasico("G1", bp=B1, Sn=750 * M, Vn=Vb13, R=0, X=0.1714))
    print(G2 := GeneradorBasico("G2", bp=B1, Sn=500 * M, Vn=Vb13, R=0, X=0.2636))
    print(G3 := GeneradorBasico("G3", bp=B1, Sn=500 * M, Vn=Vb13, R=0, X=0.2856))
    print(G4 := GeneradorBasico("G4", bp=B4, Sn=500 * M, Vn=Vb13, R=0, X=0.2513))
    print(G5 := GeneradorBasico("G5", bp=B8, Sn=500 * M, Vn=Vb13, R=0, X=0.2403))
    dt = 0.2 / 24
    print(TX1 := TransformadorTap("TX1", dt=dt, pos=0, Sn=SbTxA, Vp=Vb24, Vs=Vb13, Zev=0.0925, rel=inf, bp=B2, bs=B1))
    print(TX2 := TransformadorTap("TX2", dt=dt, pos=0, Sn=SbTxA, Vp=Vb24, Vs=Vb13, Zev=0.0925, rel=inf, bp=B2, bs=B1))
    print(TX3 := TransformadorTap("TX3", dt=dt, pos=0, Sn=SbTxB, Vp=Vb24, Vs=Vb13, Zev=0.0980, rel=inf, bp=B3, bs=B4))
    print(TX4 := TransformadorTap("TX4", dt=dt, pos=-6, Sn=SbTxB, Vp=Vb24, Vs=Vb13, Zev=0.1035, rel=inf, bp=B5, bs=B6))
    print(TX5 := TransformadorTap("TX5", dt=dt, pos=-5, Sn=SbTxB, Vp=Vb24, Vs=Vb13, Zev=0.0996, rel=inf, bp=B7, bs=B8))
    print(L1 := LineaReal("Orquidea-Alamo", Rn=0.1613, Xn=0.5197, Bn=0.94525 * u, Gn=0, bp=B2, bs=B3, L=6.5, T=4))
    print(L2 := LineaReal("Alamo-Sauce", Rn=0.1984, Xn=0.4589, Bn=0.63017 * u, Gn=0, bp=B3, bs=B5, L=13, T=3))
    print(L3 := LineaReal("Roble-Sauce", Rn=0.1984, Xn=0.4589, Bn=0.63017 * u, Gn=0, bp=B5, bs=B7, L=13.25, T=3))
    print(L4 := LineaReal("Orquidea-Roble", Rn=0.2614, Xn=0.49345, Bn=0.47263 * u, Gn=0, bp=B2, bs=B7, L=15.50, T=4))
    Y10 = G1.Y + G2.Y + G3.Y
    Y12 = TX1.Yps + TX2.Yps
    Y20 = L1.Yp0 + L4.Yp0
    Y23 = L1.Yps
    Y27 = L4.Yps
    Y30 = L1.Ys0 + L2.Yp0
    Y34 = TX3.Yps
    Y40 = G4.Y
    Y35 = L2.Yps
    Y50 = L2.Ys0 + L3.Yp0 + TX4.Yp0
    Y56 = TX4.Yps
    Y60 = TX4.Ys0
    Y57 = L3.Yps
    Y70 = L3.Ys0 + L4.Ys0 + TX5.Yp0
    Y78 = TX5.Yps
    Y80 = TX5.Ys0 + G5.Y
    Ym = np.array(
        [
            [Y10 + Y12, -Y12, 0, 0, 0, 0, 0, 0],
            [-Y12, Y20 + Y12 + Y23 + Y27, -Y23, 0, 0, 0, -Y27, 0],
            [0, -Y23, Y30 + Y34 + Y35, -Y34, -Y35, 0, 0, 0],
            [0, 0, -Y34, Y40 + Y34, 0, 0, 0, 0],
            [0, 0, 0, 0, Y50 + Y56 + Y57, -Y56, -Y57, 0],
            [0, 0, 0, 0, -Y56, Y60 + Y56, 0, 0],
            [0, -Y27, 0, 0, 0, 0, Y70 + Y27 + Y78, -Y78],
            [0, 0, 0, 0, 0, 0, -Y78, Y78 + Y80],
        ]
    )
    Ym = reduccionKron(Ym, 7)
    Ym = reduccionKron(Ym, 5)
    Ym = reduccionKron(Ym, 3)
    Ym = reduccionKron(Ym, 2)
    print()
    print("Matris de Admitancias Reducidas")
    print(Ym)
    Zm = np.linalg.inv(Ym)
    print()
    print("Matris de Impedancias Reducidas")
    print(Zm)
    C = np.array(
        [
            [1, 1, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 0, 1],
        ]
    )
    print("Matris de separación C")
    print(C)


def trabajo_final_despacho():
    G1 = GeneradorDespacho("G1", a=0.0096, b=6.4, c=400, Pmax=600, Pmin=50, Ccon=1400, Cdes=2800)
    G2 = GeneradorDespacho("G2", a=0.0080, b=8.0, c=500, Pmax=400, Pmin=50, Ccon=1600, Cdes=3200)
    G3 = GeneradorDespacho("G3", a=0.0100, b=7.9, c=600, Pmax=400, Pmin=50, Ccon=1500, Cdes=3000)
    G4 = GeneradorDespacho("G4", a=0.0110, b=7.5, c=400, Pmax=400, Pmin=50, Ccon=1450, Cdes=2900)
    G5 = GeneradorDespacho("G5", a=0.0110, b=7.5, c=400, Pmax=400, Pmin=50, Ccon=1450, Cdes=2900)
    for gen in (G1, G2, G3, G4, G5):
        print(gen)
    X1 = Combinacion("X1", G1, G2, G3, G4, G5)
    X2 = Combinacion("X2", G1, G3, G4, G5)
    X3 = Combinacion("X3", G1, G2, G4, G5)
    X4 = Combinacion("X4", G1, G2, G3, G5)
    X5 = Combinacion("X5", G1, G2, G3, G4)
    GC = GrupoCombinaciones(X1, X2, X3, X4, X5)
    print(GC)
    print(GC.reporteTransiciones())
    SL1 = CargaDespacho("SL1")
    SL2 = CargaDespacho("SL2")
    SL3 = CargaDespacho("SL3")
    B1 = BarraDespacho("B1", G1, G2, G3)
    B4 = BarraDespacho("B4", G4, carga=SL1)
    B6 = BarraDespacho("B6", carga=SL2)
    B8 = BarraDespacho("B8", G5, carga=SL3)


if __name__ == "__main__":
    trabajo_final()
    trabajo_final_despacho()
