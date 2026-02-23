import matplotlib.pyplot as plt
import numpy as np

def parameterInputs ():
    print ("Input the following parameters: \n")

    n = 1
    T0 = CR = V = k = Qin = Cycle = None

    while n <= 6:
        try:
            if n == 1:
                val = float(input("Initial temperature (T0) [K]:"))
                if val <= 0: raise ValueError("Temp must be positive.")
                T0 = val
                
            elif n == 2:
                val = float(input("Compression ratio (CR): "))
                if val <= 1: raise ValueError("CR must be greater than 1")
                CR = val

            elif n == 3:
                val = float(input("Displacement volume (V) [m^3]:"))
                if val <= 0: raise ValueError("Volume must be positive")
                V = val

            elif n == 4:
                val = float(input("Adiabatic exponent (k) :"))
                if not (1.0 < val < 1.7): raise ValueError("Not within range (normally between 1.2 and 1.4)")
                k = val

            elif n == 5:
                val = float(input("Specific heat input (Qin) [kJ/kg]:"))
                Qin = val

            elif n == 6:
                val = input("Cycle type (Otto or Diesel):").strip().lower()
                if val not in ["otto", "diesel"]: raise ValueError("Must choose between 'Otto' or 'Diesel'")
                Cycle = val

            n += 1

        except ValueError as e:
            print(f" Invalid Input: {e} try again \n")

    return T0, CR, V, k , Qin, Cycle

def cycle(T1, P1, r, Vd, k, Q, cycle_type):
    R = 0.287 # kJ/kg-k
    cv = R/(k-1)
    cp = k*cv

    v1 = (R*T1)/P1
    V1 = (v1*Vd)/(v1 - (v1/r)) #from Vd = v1-v2
    m = (P1*V1)/(R*T1)

    v2 = v1/r
    T2 = T1*(r**(k-1))
    P2 = P1*(r**k)

    if cycle_type.lower() == "otto":
        v3 = v2
        T3 = T2 + (Q/cv)
        P3 = P2 * (T3/T2)
        v4 = v1
        T4 = T3*(1/r)**(k-1)
    
    else:
        P3 = P2
        T3 = T2 + (Q/cp)
        v3 = (R*T3)/P3
        rc = v3/v2 #cutoff ratio
        v4 = v1
        T4 = T3 * (rc/r)**(k-1)

    #Work and efficiency
    W_comp = m * cv * (T1-T2) #kJ
    W_exp = m * cv * (T3-T4) if cycle_type == "otto" else m * (cp*(T3-T2)- cv*(T4-T1))
    if cycle_type == "otto":
        W_exp = m * cv *(T3-T4)
    else:
        W_exp = m *(cp*(T3-T2) - cv*(T4-T1) + (W_comp/m)) #Place holder
            