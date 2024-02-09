import cmath
from sympy import *
import numpy as np
from vpython import gcurve, color, rate, graph, label, vector, canvas


# Wprowadzanie parametrów za pomocą konsoli
U_o = float(input('Wprowadź składową stałą napięcia wyrażoną w woltach [V]: '))
f = float(input('Wprowadź częstotliwość wyrażoną w Hz (max 60 Hz) : '))
R = float(input('Wprowadź rezystancję opornika wyrażoną w Ω z przedziału 10 - 20: '))
L = float(input('Wprowadź indukcyjność cewki wyrażoną w mH z przedziału 50 - 200: '))
C = float(input('Wprowadź pojemność kondesatora wyrażoną w μF z przedziału 50 - 200: '))


# Wykres prądu i napięcia dla całego obwodu (na źródle)
canva = canvas(align='right', background=color.white)

g1 = graph(width=800, height=350, title="Wykresy prądu [*10^(-2) A] (niebieski) i napięcia [V] (czerwony) na źródle", xtitle="t [s]", ytitle="")
a = gcurve(color=color.blue, label="I(t)")
b = gcurve(color=color.red, label="U(t)")

# Wykres napięć dla kolejnych elementów obwodu RLC
g2 = graph(width=800, height=350, title="Napięcie na elementach: R - pomarańczowy, L - zielony, C - fioletowy",
           xtitle="t [s]", ytitle="U [V]")
c = gcurve(color=color.orange, label="Ur(t)")
d = gcurve(color=color.green, label="Ul(t)")
e = gcurve(color=color.purple, label="Uc(t)")


# Pulsacja i faza
omega = 2 * np.pi * f
phi = 0

# Parametry elementów układu
L_skal = L / 1000
C_skal = C / 1000000

# Wyliczenie częstotliwości rezonansowej
f_rezon = 1 / (2 * np.pi * np.sqrt(L_skal * C_skal))
print("Częstotliwość rezonansowa dla tych parametrów wynosi: " + str(f_rezon) + " Hz")

# Napięcie źródła
def U_in(t):
    return U_o * np.sin(omega * t + phi)


# Parametry czasowe do animacji
t = 0
dt = 0.001


while t < 7/f:
    rate(30)
    # Wyliczanie przesunięcia fazowego dla wykresu prądu
    fi = (omega * L_skal - 1 / (omega * C_skal)) / R
    arg = cmath.atan(fi)
    # Wyliczanie składowej stałej prądu i wzór na prąd
    Io = U_o/(np.sqrt(R**2 + (omega*L_skal - 1/(omega*C_skal))**2))
    It = Io * np.sin(omega * t - arg.real)
    # Wykres napięcia na źródle
    Ut = U_in(t)

    # Wykres napięcia na cewce
    t1 = symbols('t')
    derv_it = (Io*sin(omega * t1 - arg.real)).diff(t1)
    derv_it = lambdify(t1, derv_it)
    U_l = -(L_skal * derv_it(t))

    # Wykres napięcia na kondensatorze
    t2 = symbols('t')
    integral_it = (Io*sin(omega * t2 - arg.real)).integrate(t2)
    #integral_it = integrate(Io*sin(omega * t2 - arg.real), t2)
    integral_it = lambdify(t2, integral_it)
    U_c = -((1 / C_skal) * integral_it(t))

    # Wykres napięcia na rezystorze
    U_r = R*It
    # printing
    label(pos=vector(10, 8, 0), text=f'U_o = {U_o} V', canvas=canva)
    label(pos=vector(10, -8, 0), text=f'f = {f} Hz', canvas=canva)
    label(pos=vector(-10, 8, 0), text=f'R = {R} Ω', canvas=canva)
    label(pos=vector(-10, -8, 0), text=f'L = {L} mH', canvas=canva)
    label(pos=vector(-10, 0, 0), text=f'C = {C} μF', canvas=canva)
    label(pos=vector(10, 0, 0), text=f'f_rez = {round(f_rezon,3)} Hz', canvas=canva)
    label(pos=vector(0, 8, 0), text=f'I_o = {round(Io,3)} A', canvas=canva)
    label(pos=vector(0, 2, 0), text=f'Ur_max = {round(Io*R,3)} V ', canvas=canva)
    label(pos=vector(0, -8, 0), text=f'Ul_max = {round(omega*Io*L_skal,3)} V', canvas=canva)
    label(pos=vector(0, -2, 0), text=f'Uc_max = {round(abs((Io/omega)*(1/C_skal)),3)} V', canvas=canva)
    # Rysowanie wykresów
    a.plot((t, 100*It))
    b.plot((t, Ut))
    c.plot((t, U_r))
    d.plot((t, U_l))
    e.plot((t, U_c))
    t = t + dt
