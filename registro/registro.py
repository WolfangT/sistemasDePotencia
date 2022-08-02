# Usos y practicas

# Pip
import numpy as np

from sistemasDePotencia.redes import (
    Barra,
    Linea,
    CargaIdeal,
    GaussSiedel,
    GeneradorIdeal,
    BancoCapasitores,
    TransformadorTap,
)
from sistemasDePotencia.despacho import (
    BarraDespacho,
    GeneradorDespacho,
    tabla_costos,
    tabla_partes_iguales,
    tabla_despacho_economnico,
)

# def examen_practica():
#     Sb = 50000000
#     B4 = Barra("Barra 4", Sb=Sb, Vb=26000)
#     print(B4)
#     B8 = Barra("Barra 8", Sb=Sb, Vb=B4.Vb)
#     print(B8)
#     B9 = Barra("Barra 9", Sb=Sb, Vb=B8.Vb * 14.4 / 24)
#     print(B9)
#     B6 = Barra("Barra 6", Sb=Sb, Vb=B9.Vb * 34.5 / 13.8)
#     print(B6)
#     B7 = Barra("Barra 7", Sb=Sb, Vb=B6.Vb)
#     print(B7)
#     B2 = Barra("Barra 2", Sb=Sb, Vb=B7.Vb * 138 / 34.5)
#     print(B2)
#     B3 = Barra("Barra 3", Sb=Sb, Vb=B2.Vb)
#     print(B3)
#     B1 = Barra("Barra 1", Sb=Sb, Vb=B6.Vb * 13.8 / 34.5)
#     print(B1)
#     B5 = Barra("Barra 5", Sb=Sb, Vb=B4.Vb * 13.8 / 24)
#     print(B5)

#     L23 = Linea(
#         "Linea 2-3",
#         Sn=25000000,
#         Vn=138000,
#         Rn=8.5698,
#         Xn=51.0379,
#         Gn=0,
#         Bn=0.000315,
#         bp=B2,
#         bs=B3,
#     )
#     print(L23)

#     Tx2 = TransformadorSimple(
#         "Tx2",
#         Sn=25000000,
#         Vp=138000,
#         bp=B3,
#         Vs=24000,
#         bs=B4,
#         Zev=0.075,
#         rel=4,
#     )
#     print(Tx2)
#     Tx4 = TransformadorSimple(
#         "Tx4",
#         Sn=25000000,
#         Vp=24000,
#         bp=B8,
#         Vs=14400,
#         bs=B9,
#         Zev=0.0525,
#         rel=4,
#     )
#     print(Tx4)

#     Tx5 = TransformadorTap(
#         "Tx5",
#         Sn=25000000,
#         Vp=24000,
#         bp=B4,
#         Vs=13800,
#         bs=B5,
#         Zev=0.0475,
#         rel=4,
#         dt=0.1 / 4,
#         pos=+3,
#     )
#     print(Tx5)


# def examen_practica_2():
#     Sb = 100000000

#     B6 = Barra("Barra 6", Sb=Sb, Vb=65000)
#     B7 = Barra("Barra 7", Sb=Sb, Vb=B6.Vb)
#     B8 = Barra("Barra 8", Sb=Sb, Vb=B7.Vb * 138 / 69)
#     B9 = Barra("Barra 9", Sb=Sb, Vb=B7.Vb * 14.4 / 69)
#     B4 = Barra("Barra 4", Sb=Sb, Vb=B8.Vb)
#     B5 = Barra("Barra 5", Sb=Sb, Vb=B4.Vb * 20 / 115)
#     B3 = Barra("Barra 3", Sb=Sb, Vb=B4.Vb)
#     B2 = Barra("Barra 2", Sb=Sb, Vb=B3.Vb)
#     B1 = Barra("Barra 1", Sb=Sb, Vb=B2.Vb * 18 / 115)
#     print(B6)
#     print(B7)
#     print(B8)
#     print(B9)
#     print(B4)
#     print(B5)
#     print(B3)
#     print(B2)
#     print(B1)

#     Tx3 = TransformadorSimple(
#         "Tx3", Sn=100000000, Vp=115000, bp=B4, Vs=20000, bs=B5, Zev=0.12, rel=6.2
#     )
#     print(Tx3)
#     Tx2 = TransformadorTriple(
#         "Tx2",
#         Snp=60000000,
#         Sns=35000000,
#         Snt=25000000,
#         Vnp=138000,
#         Vns=69000,
#         Vnt=14400,
#         bp=B8,
#         bs=B7,
#         bt=B9,
#         Zept=rect(0.055, 90),
#         Zeps=rect(0.063, 90),
#         Zest=rect(0.072, 90),
#     )
#     print(Tx2)
#     Tx1 = TransformadorTriple(
#         "Tx1",
#         Snp=65000000,
#         Sns=45000000,
#         Snt=20000000,
#         Vnp=115000,
#         Vns=69000,
#         Vnt=18000,
#         bp=B2,
#         bs=B6,
#         bt=B1,
#     )
#     Tx1.pruebas_corto_circuito(Vps=3000, Ips=300, Vpt=2000, Ipt=100, Vst=1600, Ist=560)
#     print(Tx1)
#     Cs = BancoCondensadores("Cs", Qn=26000, In=500, Prel=0.013, bp=B4, bs=B3)
#     print(Cs)


# def examen_1():
#     Sb = 70000000

#     B3 = Barra("Barra 3", Sb=Sb, Vb=120000)
#     B2 = Barra("Barra 2", Sb=Sb, Vb=B3.Vb)
#     B4 = Barra("Barra 4", Sb=Sb, Vb=B3.Vb)
#     B5 = Barra("Barra 5", Sb=Sb, Vb=B4.Vb * 26 / 115)
#     B8 = Barra("Barra 8", Sb=Sb, Vb=B4.Vb)
#     B9 = Barra("Barra 9", Sb=Sb, Vb=B8.Vb * 13.8 / 115)
#     B7 = Barra("Barra 7", Sb=Sb, Vb=B8.Vb * 72 / 115)
#     B6 = Barra("Barra 6", Sb=Sb, Vb=B7.Vb)
#     B1 = Barra("Barra 1", Sb=Sb, Vb=B6.Vb * 20 / 72)
#     print(B2)
#     print(B3)
#     print(B4)
#     print(B5)
#     print(B8)
#     print(B9)
#     print(B7)
#     print(B6)
#     print(B1)

#     l1 = 170
#     L23 = Linea(
#         "Linea 2-3",
#         Sn=75000000,
#         Vn=138000,
#         Rn=l1 * 0.1,
#         Xn=l1 * 0.404,
#         Gn=0,
#         Bn=l1 * 0.292 * 10 ** -6,
#         bp=B2,
#         bs=B3,
#     )
#     print(L23)
#     L48 = Linea(
#         "Linea 4-8",
#         Sn=75000000,
#         Vn=138000,
#         Rn=l1 * 0.1,
#         Xn=l1 * 0.404,
#         Gn=0,
#         Bn=l1 * 0.292 * 10 ** -6,
#         bp=B4,
#         bs=B8,
#     )
#     print(L48)
#     l2 = 45
#     L67 = Linea(
#         "Linea 6-7",
#         Sn=45000000,
#         Vn=72000,
#         Rn=l2 * 0.18,
#         Xn=l2 * 0.5,
#         Gn=0,
#         Bn=l2 * 0.7 * 10 ** -6,
#         bp=B4,
#         bs=B8,
#     )
#     print(L67)

#     Tx3 = TransformadorTap(
#         "Tx3",
#         Sn=100000000,
#         Vp=115000,
#         bp=B4,
#         Vs=26000,
#         bs=B5,
#         Zev=0.075,
#         rel=4,
#         dt=0.0455,
#         pos=-1,
#     )
#     print(Tx3)
#     Tx2 = TransformadorTriple(
#         "Tx2",
#         Snp=75000000,
#         Sns=45000000,
#         Snt=25000000,
#         Vnp=138000,
#         Vns=72000,
#         Vnt=13800,
#         bp=B8,
#         bs=B7,
#         bt=B9,
#         Zept=rect(0.05, 90),
#         Zeps=rect(0.06, 90),
#         Zest=rect(0.07, 90),
#     )
#     print(Tx2)
#     Tx1 = TransformadorTriple(
#         "Tx1",
#         Snp=75000000,
#         Sns=45000000,
#         Snt=25000000,
#         Vnp=138000,
#         Vns=72000,
#         Vnt=20000,
#         bp=B2,
#         bs=B6,
#         bt=B1,
#     )
#     Tx1.pruebas_corto_circuito(Vps=2800, Ips=280, Vpt=1800, Ipt=90, Vst=1440, Ist=504)
#     print(Tx1)

#     print(V2_fic := B6.Vb / 72 * 138)
#     print(t := B2.Vb / V2_fic)
#     print(tap_fic(Tx1.Zp, t))

#     G1 = Generador("G1", Sn=65000000, Vn=20000, Zev=rect(1.35, 90), b=B1)
#     G2 = Generador("G2", Sn=100000000, Vn=26000, Zev=rect(1.65, 90), b=B5)
#     print(G1)
#     print(G2)

#     Cs = BancoCondensadores("Cs", Qn=0.007 * 10 ** 6, In=256, Prel=0.0165, bp=B3, bs=B4)
#     print(Cs)

#     # Cs = BancoCondensadores("Cs", Qn=26000, In=500, Prel=0.013, bp=B4, bs=B3)
#     # print(Cs)


# def ejercicio_viernes():
#     Sb = 100 * 10 ** 6
#     Vb = 230000
#     B5 = Barra("Tablaso", Sb=Sb, Vb=Vb)
#     B4 = Barra("Morochas", Sb=Sb, Vb=Vb)
#     print(B5)
#     print(B4)
#     L45 = Linea(
#         "L45",
#         L=67,
#         Rn=0.0702,
#         Xn=0.465000j,
#         Gn=0,
#         Bn=3.473j * 10 ** -6,
#         In=850,
#         bp=B5,
#         bs=B4,
#     )
#     print(L45)
#     print(V4 := V("V4", rect(0.9662, -16.4), B5))
#     print(V5 := V("V5", rect(0.98636, -11.8), B5))
#     print(I1 := I("I1", V5 * L45.B / 2, B5))
#     print(I2 := I("I2", V4 * L45.B / 2, B5))
#     print(I3 := I("I3", (V5 - V4) / (L45.R + L45.X), B5))

#     print(Sent := S("S_ent", V5 * (I1 + I3).conjugate(), B5))
#     print(Ssal := S("S_sal", V4 * (I3 - I2).conjugate(), B5))
#     print(S1 := S("S1", V5 * I1.conjugate(), B5))
#     print(S2 := S("S2", V4 * I2.conjugate(), B5))

#     print(Sent.en_real())
#     print(Ssal.en_real())
#     print(S1.en_real())
#     print(S2.en_real())


# def clase_19_03():
#     Sb = 100 * 10 ** 6
#     Vb = 230000
#     B1 = Barra("B1", Sb=Sb, Vb=Vb)
#     B2 = Barra("B2", Sb=Sb, Vb=Vb)
#     B3 = Barra("B3", Sb=Sb, Vb=Vb)
#     B4 = Barra("B4", Sb=Sb, Vb=Vb)
#     B5 = Barra("B5", Sb=Sb, Vb=Vb * 115 / 230)
#     print(B1)
#     print(B2)
#     print(B3)
#     print(B4)
#     print(B5)
#     L12 = Linea("L12", R=0.01008, X=0.05040, G=0, B=0.05125 * 2, bp=B1, bs=B2)
#     print(L12)
#     # print(L12.flujo_potencias())
#     L13 = Linea("L13", R=0.00744, X=0.03720, G=0, B=0.03875 * 2, bp=B1, bs=B3)
#     print(L13)
#     # print(L13.flujo_potencias())
#     L24 = Linea("L24", R=0.00744, X=0.03720, G=0, B=0.03875 * 2, bp=B2, bs=B4)
#     print(L24)
#     # print(L24.flujo_potencias())
#     L34 = Linea("L34", R=0.01272, X=0.06360, G=0, B=0.06375 * 2, bp=B3, bs=B4)
#     print(L34)
#     # print(L34.flujo_potencias())
#     TX1 = TransformadorTap(
#         "TX1",
#         dt=(1 - 0.96385542),
#         pos=-1,
#         Sn=300 * 10 ** 6,
#         Vp=230000,
#         Vs=115000,
#         Zev=0.06,
#         rel=float("inf"),
#         bp=B3,
#         bs=B5,
#     )
#     print(TX1)
#     # print(TX1.flujo_potencias())
#     print(Ic := I("Ic", (-52.3j * 10 ** 6 / B5.Vb).conjugate()))
#     print(Yc := Y("Yc", (Ic / B5.Vb)))
#     print(Ycpu := Y("Yc", Yc / B5.Yb, B5))
#     Ym = np.matrix(
#         [
#             [L12.Yp0 + L12.Yps + L13.Yp0 + L13.Yps, -L12.Yps, -L13.Yps, 0, 0],
#             [-L12.Yps, L12.Ys0 + L12.Yps + L24.Yp0 + L24.Yps, 0, -L24.Yps, 0],
#             [
#                 -L13.Yps,
#                 0,
#                 L13.Yp0 + L13.Ys0 + L34.Yp0 + L34.Yps + TX1.Yps + TX1.Yp0,
#                 -L34.Yps,
#                 -TX1.Yps,
#             ],
#             [0, -L24.Yps, -L34.Yps, L24.Yps + L24.Ys0 + L34.Yps + L34.Ys0, 0],
#             [0, 0, -TX1.Yps, 0, TX1.Yps + TX1.Ys0 + Ycpu],
#         ]
#     )
#     print(Ym)


def wolfang():
    print("Datos Wolfang")
    Sb = 100 * 10 ** 6
    Vb1 = 230000
    Vb5 = 161000
    alpha = 1.6
    print(B1 := Barra("B1", Sb=Sb, Vb=Vb1))
    B2 = Barra("B2", Sb=Sb, Vb=Vb1)
    B3 = Barra("B3", Sb=Sb, Vb=Vb1)
    B4 = Barra("B4", Sb=Sb, Vb=Vb1)
    print(B5 := Barra("B5", Sb=Sb, Vb=Vb5))
    print(L12 := Linea("L12", R=0.01008, X=0.05040, G=0, B=0.05125 * 2, bp=B1, bs=B2))
    print(L13 := Linea("L13", R=0.00744, X=0.03720, G=0, B=0.03875 * 2, bp=B1, bs=B3))
    print(L24 := Linea("L24", R=0.00744, X=0.03720, G=0, B=0.03875 * 2, bp=B2, bs=B4))
    print(L34 := Linea("L34", R=0.01272, X=0.06360, G=0, B=0.06375 * 2, bp=B3, bs=B4))
    print(L12.en_real())
    print(L13.en_real())
    print(L24.en_real())
    print(L34.en_real())
    print(
        TX1 := TransformadorTap(
            "TX1",
            dt=0.2 / 32,
            pos=+10,
            Sn=310 * 10 ** 6,
            Vp=Vb5,
            Vs=Vb1,
            Zev=0.0625,
            rel=3.487,
            bp=B5,
            bs=B3,
        )
    )
    print(
        Cx := BancoCapasitores(
            "CX1",
            bp=B5,
            Sn=47.5j * 10 ** 6 / B5.Sb,
            Vn=TX1.Vp / B5.Vb,
            pj=0.034,
        )
    )
    print(
        G2 := GeneradorIdeal(
            "G2",
            bp=B4,
            Pn=326 * 10 ** 6,
            Vn=234875.7,
            Qnmin=-6.3 * 10 ** 6,
            Qnmax=233 * 10 ** 6,
        )
    )
    print(CargaIdeal("C1", bp=B1, Pn=50 * 10 ** 6, Qn=30.99 * 10 ** 6))
    print(C2 := CargaIdeal("C2", bp=B2, Pn=170 * 10 ** 6, Qn=105.35 * 10 ** 6))
    print(C3 := CargaIdeal("C3", bp=B3, Pn=0, Qn=0))
    print(C4 := CargaIdeal("C4", bp=B4, Pn=80 * 10 ** 6, Qn=49.58 * 10 ** 6))
    print(C5 := CargaIdeal("C5", bp=B5, Pn=200 * 10 ** 6, Qn=123.94 * 10 ** 6))
    print("Matris Admitancias")
    print(
        Ym := np.array(
            [
                [L12.Yp0 + L12.Yps + L13.Yp0 + L13.Yps, -L12.Yps, -L13.Yps, 0, 0],
                [-L12.Yps, L12.Ys0 + L12.Yps + L24.Yp0 + L24.Yps, 0, -L24.Yps, 0],
                [
                    -L13.Yps,
                    0,
                    L13.Yps + L13.Ys0 + L34.Yps + L34.Yp0 + TX1.Yps + TX1.Ys0,
                    -L34.Yps,
                    -TX1.Yps,
                ],
                [0, -L24.Yps, -L34.Yps, L24.Yps + L24.Ys0 + L34.Yps + L34.Ys0, 0],
                [0, 0, -TX1.Yps, 0, TX1.Yps + TX1.Yp0 + Cx.Yn],
            ]
        )
    )
    print()
    print(
        GS := GaussSiedel(
            [B1, B2, B3, B4, B5],
            [1 + 0j, None, None, G2.Vn, None],
            [None, -C2.S.real, -C3.S.real, G2.Pn - C4.S.real, -C5.S.real],
            [None, -C2.S.imag * 1j, -C3.S.imag * 1j, None, -C5.S.imag * 1j],
            alpha,
            Ym,
        )
    )
    print(GS.iteracion())
    print(GS.iteracion())


# Potencia 2
def evaluacion_1_ejercicio_1():
    barra = BarraDespacho(
        "B1",
        GeneradorDespacho("1", a=0.012, b=9),
        GeneradorDespacho("2", a=0.0096, b=6),
        GeneradorDespacho("3", a=0.008, b=8),
        GeneradorDespacho("4", a=0.0068, b=10),
    )
    print(barra)
    data = barra.analisis(80)
    print("Partes Iguales")
    print(tabla_partes_iguales(data))
    print()
    print("Despacho Economnico")
    print(tabla_despacho_economnico(data))
    print()
    print("Costos")
    print(tabla_costos(data))
    print()


if __name__ == "__main__":
    wolfang()
