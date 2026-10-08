# Three MATLAB Models: Code and LaTeX Formulae

---

## 1. Solar PV Module (Single-Diode Model)

### Code

```matlab
clc

k = 1.38065e-23; q = 1.602e-19;
Isc = 8.21; Voc = 32.9; Ki = 0.0032;
Ns = 54; T = 298; Tn = 303;
G = 1000; Gn = 1000;
a = 2; Eg = 1.2;
Rs = 0.221; Rp = 415.405;

Vtn = Ns*k*Tn/q;
Vt = Ns*k*T/q;
I0n = Isc/(exp(Voc/(a*Vtn)) - 1);
I0 = I0n*(Tn/T)^3*exp(q*Eg/(a*k)*(1/Tn - 1/T));
Iph = G/Gn*(Isc + Ki*(T - Tn));

V = Voc:-0.1:0;
I = zeros(size(V));
Iprev = 0;

for n = 1:numel(V)
    I(n) = max(Iprev, 0);
    x = V(n) + Iprev*Rs;
    Iprev = Iph - I0*(exp(x/(a*Vt)) - 1) - x/Rp;
end

plot(V, I, 'r', 'LineWidth', 2.5)
xlabel('Voltage (V)')
ylabel('Current (A)')

figure
plot(V, V.*I, 'k', 'LineWidth', 2.5)
xlabel('Voltage (V)')
ylabel('Power (W)')
```

### Formulae

Thermal voltage (actual and nominal):

$$V_t = \frac{N_s\,k\,T}{q}, \qquad V_{t,n} = \frac{N_s\,k\,T_n}{q}$$

Nominal saturation current:

$$I_{0,n} = \frac{I_{sc}}{\exp\left(\dfrac{V_{oc}}{a\,V_{t,n}}\right) - 1}$$

Saturation current at temperature $T$:

$$I_0 = I_{0,n}\left(\frac{T_n}{T}\right)^{3}\exp\left[\frac{q\,E_g}{a\,k}\left(\frac{1}{T_n} - \frac{1}{T}\right)\right]$$

Photocurrent:

$$I_{ph} = \frac{G}{G_n}\left[I_{sc} + K_i\,(T - T_n)\right]$$

Terminal current (implicit):

$$I = I_{ph} - I_0\left[\exp\left(\frac{V + I\,R_s}{a\,V_t}\right) - 1\right] - \frac{V + I\,R_s}{R_p}$$

Power:

$$P = V\,I$$

Symbols: $k$ Boltzmann constant, $q$ electron charge, $N_s$ series cells, $a$ diode ideality factor, $E_g$ silicon band gap, $R_s$ series resistance, $R_p$ parallel resistance, $I_{sc}$ short-circuit current, $V_{oc}$ open-circuit voltage, $K_i$ current temperature coefficient, $G$ irradiance, $T$ temperature in kelvin.

---

## 2. PEM Fuel Cell Stack (Polarization Curve)

### Code

```matlab
clear; clc; close all

T = 310; pH2 = 1; pO2 = 1;
A = 69.7; N = 24;
z1 = -0.475; z3 = 7.6e-5; z4 = -1e-4;
Rc = 0.00019; Rm = 0.2;
B = 0.0171; Jmax = 1600;

I = 0.1:0.1:9.9;
J = I/A;

E = 1.229 - 0.85e-3*(T - 298.15) + 1.31e-5*T*(log(pH2) + 0.5*log(pO2));
cO2 = pO2/(5.08e6*exp(-498/T));
z2 = 0.00286 + 0.0002*log(A) + 4.3e-5*log(cO2);

Vact = -(z1 + z2*T + z3*T*log(cO2) + z4*T*log(I));
Vohm = I*(Rm + Rc);
Vcon = -B*log(1 - J/Jmax);
V = N*(E - Vact - Vohm - Vcon);

yyaxis left
plot(I, V)
ylabel('Voltage (V)')
yyaxis right
plot(I, V.*I)
ylabel('Power (W)')
xlabel('Current (A)')
```

### Formulae

Nernst voltage:

$$E = 1.229 - 0.85\times10^{-3}\,(T - 298.15) + 1.31\times10^{-5}\,T\left[\ln p_{H_2} + 0.5\ln p_{O_2}\right]$$

Oxygen concentration:

$$c_{O_2} = \frac{p_{O_2}}{5.08\times10^{6}\,\exp\left(-498/T\right)}$$

Activation coefficient:

$$z_2 = 0.00286 + 0.0002\ln A + 4.3\times10^{-5}\ln c_{O_2}$$

Activation loss:

$$V_{act} = -\left[z_1 + z_2\,T + z_3\,T\ln c_{O_2} + z_4\,T\ln I\right]$$

Ohmic loss:

$$V_{ohm} = I\,(R_m + R_c)$$

Current density:

$$J = \frac{I}{A}$$

Concentration loss:

$$V_{con} = -B\ln\left(1 - \frac{J}{J_{max}}\right)$$

Stack voltage:

$$V = N\left(E - V_{act} - V_{ohm} - V_{con}\right)$$

Stack power:

$$P = V\,I$$

Symbols: $T$ stack temperature in kelvin, $p_{H_2}$ and $p_{O_2}$ partial pressures in atm, $A$ cell active area, $N$ number of cells, $z_1$ to $z_4$ empirical activation coefficients, $R_m$ membrane resistance, $R_c$ contact resistance, $B$ concentration-loss constant, $J_{max}$ limiting current density.

---

## 3. Lead-Acid Battery (Charging Simulation)

### Code

```matlab
clear; clc; close all

I = 5;
K = 0.8; D = 1e-5;
Cap = 936; ns = 6;
SOC0 = 0.2;
t = 0:0.1:7;

SOC = zeros(size(t));
Vbat = zeros(size(t));
s = SOC0;

for n = 1:numel(t)
    if I <= 0
        Voc = (1.926 + 0.124*s)*ns;
        R = (0.19 + 0.1037/(s - 0.14))*ns/Cap;
    else
        Voc = (2 + 0.148*s)*ns;
        R = (0.758 + 0.1309/(1.06 - s))*ns/Cap;
    end
    s = SOC0 + (K*Voc*I - D*s*Cap)*t(n)/Cap;
    SOC(n) = s;
    Vbat(n) = Voc + I*R;
end

Vbat
SOC
plot(t, Vbat)
figure
plot(t, SOC)
```

### Formulae

Open-circuit voltage:

$$V_{oc} = \begin{cases} n_s\,(2 + 0.148\,SOC), & I > 0 \\[4pt] n_s\,(1.926 + 0.124\,SOC), & I \le 0 \end{cases}$$

Internal resistance:

$$R = \frac{n_s}{C}\begin{cases} 0.758 + \dfrac{0.1309}{1.06 - SOC}, & I > 0 \\[8pt] 0.19 + \dfrac{0.1037}{SOC - 0.14}, & I \le 0 \end{cases}$$

Net power:

$$P_{net} = K\,V_{oc}\,I - D\,SOC\,C$$

State of charge:

$$SOC = SOC_0 + \frac{P_{net}\,t}{C}$$

Terminal voltage:

$$V_{bat} = V_{oc} + I\,R$$

Symbols: $I$ battery current (positive charges), $K$ charging efficiency, $D$ self-discharge rate, $C$ capacity (`Cap` in the code), $n_s$ cells in series, $SOC_0$ initial state of charge, $t$ time.
