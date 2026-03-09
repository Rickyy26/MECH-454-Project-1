import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

R_gas = 0.287 #kJ/kg-K
L_rod = 0.15 #Connection rod length assumption [m]
a_crank = 0.045 # crank radius assumption [m]

def full_cycle(T1, P1, Vd, r, k, Q, N, cycle_type):
    cv = R_gas/(k-1)
    cp = k * cv

    #initial state
    V1_total = Vd/(1-(1/r))
    m = (P1*V1_total)/(R_gas*T1)

    #compression (1->2)
    theta_comp = np.linspace(np.pi, 0, 100)
    T_comp = odeint(engine_ode, T1, theta_comp, args=(m, Vd, r, k, R_gas, L_rod, a_crank, 0))
    T2 = T_comp[-1][0]

    #heat addition (2->3)
    if cycle_type == 'otto':
        T3 = T2 + (Q/cv)
        V3 = V1_total/r

    else:
        T3 = T2 + (Q/cp)
        V3 = (m*R_gas*T3)/(P1 * (r**k)) #constant P @ v3

    #expensaion (3->4)
    v_ratio = V3/(V1_total/r)
    r_exp = r/v_ratio
    theta_exp = np.linspace(0, np.pi, 100)
    T_exp = odeint(engine_ode, T3, theta_exp, args=(m, Vd, r_exp, k, R_gas, L_rod, a_crank, 0))
    T4 = T_exp[-1][0]

    #parameters
    Qin_total = m * Q
    Qout_total = m * cv * (T4-T1)
    W_net = Qin_total - Qout_total
    power = (W_net * N)/120 #kW for a 4 stroke
    torque = (W_net * 1000)/(4 * np.pi) #Nm
    eta = W_net/Qin_total
    mep = W_net/Vd
    spec_power = power/(np.pi * (0.124**2)/4 * 0.124 * 1000) #kW/L placeholder

    return{
        "power": power,
        "torque": torque,
        "efficiency": eta,
        "mep": mep,
    }

def crank(theta, Vd, r, L, a):
    R=L/a
    V_c = Vd/(r-1)

    term1 = 0.5*(r-1)
    term2 = R + 1 - np.cos(theta) - np.sqrt(R**2 - np.sin(theta)**2)
    V = V_c * (1 + term1 * term2)

    dV_dtheta = V_c * term1 * (np.sin(theta) + (np.sin(theta)*np.cos(theta))/np.sqrt(R**2 - np.sin(theta)**2))
    return V, dV_dtheta

def engine_ode(T, theta, m, Vd, r, k, R_gas, L, a, dQ_dtheta):
    cv = R_gas/(k-1)
    V, dV_dtheta = crank(theta, Vd, r, L, a)
    P = (m * R_gas * T)/V

    dT_dtheta = (1/(m * cv)) * (dQ_dtheta - P * dV_dtheta)
    return dT_dtheta

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

def isentropic(T_start, P_start, v_start, v_end, Vd, r, k, npts):
    R = 0.287 #kJ/kg-k
    v_vec = np.linspace(v_start, v_end, npts)
    P_vec = P_start * (v_start/v_vec)**k
    T_vec = T_start * (v_start/v_vec)**(k-1)

    #work integral
    P_end = P_vec[-1]
    specific_work = (P_end * v_end - P_start * v_start)/(1-k)
    V1_total = Vd/(1-(1/r))
    mass = (P_start*V1_total)/(R*T_start) if v_start > v_end else (P_vec[0]*V1_total)/(R*T_vec[0])
    m = (P_start*(v_start*mass/v_start))/(R*T_start) #placeholder

    return T_vec, P_vec, v_vec, specific_work * mass

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
        P4 = P3*(v3/v4)**k
    
    else:
        P3 = P2
        T3 = T2 + (Q/cp)
        v3 = (R*T3)/P3
        rc = v3/v2 #cutoff ratio
        v4 = v1
        T4 = T3 * (rc/r)**(k-1)
        P4 = P3*(v3/v4)**k

    #Work and efficiency calculated
    W_comp = m * cv * (T1-T2) #kJ
    Q_out = m * cv * (T4-T1) #kJ
    W_net = (m * Q) - Q_out #kJ
    W_exp = W_net - W_comp #kJ
    eta = W_net/(m*Q)
    torque = (W_net*1000)/(4*np.pi) #Nm

    states = {
        "st1":{"T": T1, "P": P1, "v": v1},
        "st2":{"T": T2, "P": P2, "v": v2},
        "st3":{"T": T3, "P": P3, "v": v3},
        "st4":{"T": T4, "P": P4, "v": v4}
    }
    return abs(W_comp), W_exp, W_net, eta, torque, states


#       ---MAIN CODE---

#Q4

#engine specs given
Vd_liters = 1.5
Vd_m3 = Vd_liters/1000
Q_input = 1800
T1_input = 15+273.15
P1_input = 100
k_input = 1.4

rpm_range = np.linspace(1000, 6000, 50)

#simulations
otto_results = [full_cycle(T1_input, P1_input, Vd_m3, 10, k_input, Q_input, n, 'otto') for n in rpm_range]
diesel_results = [full_cycle(T1_input, P1_input, Vd_m3, 20, k_input, Q_input, n, 'diesel') for n in rpm_range]

#plotting
fig, axs = plt.subplots(2,2,figsize=(12,10))

#power plot
axs[0,0].plot(rpm_range, [d['power'] for d in otto_results], 'r', label='Otto')
axs[0,0].plot(rpm_range, [d['power'] for d in diesel_results], 'b', label='Diesel')
axs[0,0].set_title("Power vs RPM"); axs[0,0].set_ylabel("Power [kW]"); axs[0,0].legend()

#torque plt
axs[0,1].plot(rpm_range, [d['torque'] for d in otto_results], 'r')
axs[0,1].plot(rpm_range, [d['torque'] for d in diesel_results], 'b')
axs[0,1].set_title("Torque vs RPM"); axs[0,1].set_ylabel("Torque [Nm]")

#efficiency plot
axs[1,0].plot(rpm_range, [d['efficiency'] for d in otto_results], 'r')
axs[1,0].plot(rpm_range, [d['efficiency'] for d in diesel_results], 'b')
axs[1,0].set_title("Thermal efficiency vs RPM"); axs[1,0].set_ylabel("Efficiency")
axs[1,0].set_ylim(0,1)

#mep plot
axs[1,1].plot(rpm_range, [d['mep'] for d in otto_results], 'r')
axs[1,1].plot(rpm_range, [d['mep'] for d in diesel_results], 'b')
axs[1,1].set_title("MEP vs RPM"); axs[1,1].set_ylabel("Mep [kPa]")

for ax in axs.flat: ax.set_xlabel("Engine speed [RPM]"); ax.grid(True)
plt.tight_layout(); plt.show()

#Q1 constants
T1_const = 15 + 273 #K
P1_const = 90 #kPa
Vd_const = 1.5/1000 #m^3
Q_const = 1800 #kJ/kg
k_const = 1.4

T1_in, r_in, Vd_in, k_in, Q_in, cycle_in = parameterInputs()
Wc_q1, We_q1, Wn_q1, et_q1, trq_q1, sts = cycle(T1_in, 90, r_in, Vd_in, k_in, Q_in, cycle_in)

#Q3
R_gas = 0.287
V1_total = Vd_in/(1-(1/r_in))
m_total = (sts['st1']['P'] * V1_total)/(R_gas*sts['st1']['T'])
theta_comp = np.linspace(np.pi, 0, 100)     #starting at BDC and going to TDC (0)
T_num_comp = odeint(engine_ode, sts['st1']['T'], theta_comp, args=(m_total, Vd_in, r_in, k_in, R_gas, L_rod, a_crank, 0))
T2_ode = T_num_comp[-1][0]
T2_analytical = sts['st2']['T']
error_q3 = abs(T2_ode - T2_analytical)/T2_analytical*100

print(f"\n Question 3")
print(f"Analytical T2 (Q1): {T2_analytical:.4f} K")
print(f"Numberical T2 (Q3): {T2_ode:.4f} K")
print(f"Difference: {error_q3:.2e} %")

#isentropic path for q2
#Compression from 1 -> 2
T_c, P_c, v_c, W_c_int = isentropic(sts['st1']['T'], sts['st1']['P'], sts['st1']['v'], sts['st2']['v'], Vd_in, r_in, k_in, 100)
#expansion from 3 -> 4
T_e, P_e, v_e, W_e_int = isentropic(sts['st3']['T'], sts['st3']['P'], sts['st3']['v'], sts['st4']['v'], Vd_in, r_in, k_in, 100)

cases = [(288.15, 90, 8), (288.15, 100, 8), (300.15, 90, 8), (288.15, 90, 10), (288.15, 90, 12), (310.15, 95, 9)]
fig_q2, (ax_pv, ax_tv) = plt.subplots(1,2,figsize=(15,6))
print(f"\n{'Case':<5} | {'Int W [kJ]':<12} | {'State W [kJ]':<12} | {'Error %':<8}")
print("-" * 50)

for i, (t1, p1, r) in enumerate(cases): #case comparison loop
    v1_c = (0.287*t1)/p1
    v2_c = v1_c/r
    T_v, P_v, v_v, W_int = isentropic(t1, p1, v1_c, v2_c, Vd_in, r, k_in, 100)

    t2_s = t1*(r**(k_in-1))
    m_s = (p1* (Vd_in/(1-(1/r))))/(0.287*t1)
    W_state = m_s * (0.287/(k_in-1))*(t1-t2_s)

    label = f"Case {i+1}: r={r}"
    ax_pv.plot(v_v, P_v, label=label)
    ax_tv.plot(v_v, T_v, label=label)
    print(f"{i+1:<5} | {W_int:<12.4f} | {W_state:<12.4f} | {abs((W_int-W_state)/W_state)*100:<8.2e}")


ax_pv.set_title("Q2 plot 1: P vs v (6 cases)"); ax_pv.set_ylabel("P [kPa]"); ax_pv.legend(); ax_pv.grid(True)
ax_tv.set_title("Q2 plot 2: T vs v (6 cases)"); ax_tv.set_ylabel("T [K]"); ax_tv.legend(); ax_tv.grid(True)

r_range = np.linspace(1.1, 20, 50)
otto_results = [cycle(T1_const, P1_const, r, Vd_const, k_const, Q_const, "otto") for r in r_range]
diesel_results = [cycle(T1_const, P1_const, r, Vd_const, k_const, Q_const, "diesel") for r in r_range]

fig1, axs = plt.subplots(2, 2, figsize=(12, 10))
axs[0,0].plot(r_range, [d[2] for d in otto_results], 'r', label="Otto"); axs[0,0].plot(r_range, [d[2] for d in diesel_results], 'b', label="Diesel")
axs[0,0].set_title("Net Work [kJ]"); axs[0,0].legend(); axs[0,0].grid(True)
axs[0,1].plot(r_range, [d[3] for d in otto_results], 'r'); axs[0,1].plot(r_range, [d[3] for d in diesel_results], 'b')
axs[0,1].set_title("Efficiency"); axs[0,1].grid(True)
axs[1,0].plot(r_range, [d[4] for d in otto_results], 'r'); axs[1,0].plot(r_range, [d[4] for d in diesel_results], 'b')
axs[1,0].set_title("Torque [Nm]"); axs[1,0].grid(True)
axs[1,1].plot(r_range, [d[5]['st2']['T'] for d in otto_results], 'r', label="Otto"); axs[1,1].plot(r_range, [d[5]['st2']['T'] for d in diesel_results], 'b--', label="Diesel")
axs[1,1].set_title("T2 [K]"); axs[1,1].legend(); axs[1,1].grid(True)
plt.tight_layout()

#plotting Q2 (p-v diagram)
plt.figure(figsize=(10,6))
plt.plot(v_c, P_c, 'b', label='Compression (Isentropic)')
plt.plot(v_e, P_e, 'r', label='Expansion (Isentropic)')

#p-v
if cycle_in == 'otto':
    plt.plot([sts['st2']['v'], sts['st3']['v']], [sts['st2']['P'], sts['st3']['P']], 'g', label='Combustion (v=const)')
else:
    plt.plot([sts['st2']['v'], sts['st3']['v']], [sts['st2']['P'], sts['st3']['P']], 'g', label='Combustion (P=const)')

plt.plot([sts['st4']['v'], sts['st1']['v']], [sts['st4']['P'], sts['st1']['P']], 'k', label='Exhaust (v=const)')
plt.title(f"P-v Diagram: {cycle_in.capitalize()} Cycle")
plt.xlabel("Specific Volume [m^3/kg]"); plt.ylabel("Pressure [kPa]"); plt.legend(); plt.grid(True)
plt.show()

#plotting Q2 (t-v diagram)
plt.figure(figsize=(10,6))
plt.plot(v_c, T_c, 'b', label='Compression (Isentropic)')
plt.plot(v_e, T_e, 'r', label='Expansion (Isentropic)')

#t-v
if cycle_in == 'otto':
    plt.plot([sts['st2']['v'], sts['st3']['v']], [sts['st2']['T'], sts['st3']['T']], 'g', label='Combustion (v=const)')
else:
    plt.plot([sts['st2']['v'], sts['st3']['v']], [sts['st2']['T'], sts['st3']['T']], 'g', label='Combustion (P=const)')

plt.plot([sts['st4']['v'], sts['st1']['v']], [sts['st4']['T'], sts['st1']['T']], 'k', label='Exhaust (v=const)')
plt.title(f"Tt-v Diagram: {cycle_in.capitalize()} Cycle")
plt.xlabel("Specific Volume [m^3/kg]"); plt.ylabel("Temperature [K]"); plt.legend(); plt.grid(True)
plt.show()

#printing table of main points
print(f"\n {'State':<10} | {'Temp [K]':<10} | {'Pres [kPa]':<12} | {'Vol [m^3/kg]':<12}")
print("-" * 55)

for key in ['st1', 'st2', 'st3', 'st4']:
    T = sts[key]['T']
    P = sts[key]['P']
    v = sts[key]['v']
    print(f"{key:<10} | {T:<10.2f} | {P:<12.2f} | {v:<12.4f}")

#result printing
print(f"\n--- Results for {cycle_in.upper()} ---")
print(f"Efficiency: {et_q1:.2%}")
print(f"Net Work: {Wn_q1:.2f} kJ")
print(f"Compression Work: {Wc_q1:.2f} kJ")
print(f"Expansion Work: {We_q1:.2f} kJ")
print(f"Torque: {trq_q1:.2f} Nm")

r_range = np.linspace(1.1, 20, 50)
otto_results = [cycle(T1_const, P1_const, r, Vd_const, k_const, Q_const, "otto") for r in r_range]
diesel_results = [cycle(T1_const, P1_const, r, Vd_const, k_const, Q_const, "diesel") for r in r_range]

fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs[0,0].plot(r_range, [d[2] for d in otto_results], 'r', label="Otto")
axs[0,0].plot(r_range, [d[2] for d in diesel_results], 'b', label="Diesel")
axs[0,0].set_title("Net Work [kJ]"); axs[0,0].legend()

axs[0,1].plot(r_range, [d[3] for d in otto_results], 'r')
axs[0,1].plot(r_range, [d[3] for d in diesel_results], 'b')
axs[0,1].set_title("Efficiency")

axs[1,0].plot(r_range, [d[4] for d in otto_results], 'r')
axs[1,0].plot(r_range, [d[4] for d in diesel_results], 'b')
axs[1,0].set_title("Torque [Nm]")

axs[1,1].plot(r_range, [d[5]['st2']['T'] for d in otto_results], 'r')
axs[1,1].plot(r_range, [d[5]['st2']['T'] for d in diesel_results], 'b')
axs[1,1].set_title("T2 before combustion [K]") #Values are overlapping therefore graph only shows one

for ax in axs.flat: ax.set_xlabel("Compression Ratio (r)"); ax.grid(True)
plt.tight_layout(); plt.show()