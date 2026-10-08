# Lab Programs: Code and LaTeX Formulae
 
---
 
## 1. Exp7: IV PV (Simple Approximation)
 
### Code
 
```matlab
clc;
clear;
Voc = 30; % Open circuit voltage (V)32.9
Isc = 8; % Short circuit current (A)8.21
V = 0:0.1:Voc;
% Simple approximate I-V curve
I = Isc * (1 - (V/Voc).^2);
% Power
P = V .* I;
% I-V curve
figure;
plot(V,I,'LineWidth',2);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('PV I-V Characteristics');
grid on;
% P-V curve
figure;
plot(V,P,'LineWidth',2);
xlabel('Voltage (V)');
ylabel('Power (W)');
title('PV P-V Characteristics');
grid on;
```
 
### Formulae
 
Current:
 
$$I = I_{sc}\left[1 - \left(\frac{V}{V_{oc}}\right)^{2}\right]$$
 
Power:
 
$$P = V\,I$$
 
Symbols: $I_{sc}$ short-circuit current, $V_{oc}$ open-circuit voltage.
 
---
 
## 2. Exp7: IV PV (With Knee)
 
### Code
 
```matlab
clc;
clear;
close all;
Voc = 32.9; % Open circuit voltage (V)32.9
Isc = 8.21; % Short circuit current (A)8.21
V = 0:0.1:Voc;
% PV I-V characteristic
Vknee = 24; % Knee voltage
I = Isc * (1 - (V/Vknee).^10);
% Set current to zero after knee
I(V > Voc) = 0;
I(I < 0) = 0;
% Power
P = V .* I;
% I-V characteristic
figure;
plot(V,I,'LineWidth',2);xlabel('Voltage (V)');
ylabel('Current (A)');
title('PV I-V Characteristics');
grid on;
% P-V characteristic
figure;
plot(V,P,'LineWidth',2);
xlabel('Voltage (V)');
ylabel('Power (W)');
title('PV P-V Characteristics');
grid on;
```
 
### Formulae
 
Current with clamping at zero:
 
$$I = \max\left\{0,\; I_{sc}\left[1 - \left(\frac{V}{V_{knee}}\right)^{10}\right]\right\}$$
 
Power:
 
$$P = V\,I$$
 
Symbols: $I_{sc}$ short-circuit current, $V_{knee}$ knee voltage.
 
---
 
## 3. Exp6: Energy Storage System (Battery)
 
### Code
 
```matlab
clc;
clear;
% Battery parameters
Capacity = 936; % Battery capacity
SOC = 0.2; % Initial SOC = 20%
I = 5; % Charging current (A)
Vnom = 12; % Nominal battery voltage (V)
% Time
t = 0:0.1:80;
% Calculate SOC
SOC = SOC + (I * t) / Capacity;
% Limit SOC to 100%
SOC(SOC > 1) = 1;
% Simple battery voltage model
Vbat = Vnom + 2*SOC;
% Plot battery voltage
figure;
plot(t,Vbat,'LineWidth',2);
xlabel('Time (s)');
ylabel('Battery Voltage (V)');
title('Battery Voltage vs Time');
grid on;
% Plot SOC
figure;
plot(t,SOC,'LineWidth',2);
xlabel('Time (s)');
ylabel('State of Charge');
title('SOC vs Time');
grid on;
```
 
### Formulae
 
State of charge, limited to 100%:
 
$$SOC = \min\left\{1,\; SOC_0 + \frac{I\,t}{C}\right\}$$
 
Battery voltage:
 
$$V_{bat} = V_{nom} + 2\,SOC$$
 
Symbols: $SOC_0$ initial state of charge, $I$ charging current, $t$ time, $C$ battery capacity (`Capacity` in the code), $V_{nom}$ nominal battery voltage.
 
---
 
## 4. Fuel Cell
 
### Code
 
```matlab
clc;
I = 0.1:0.1:10;
T = 300;
A = 70;
z1 = -0.5;
z2=0.005;
z3 = 7e-7;
z4 = -7e-7;
co2=7e-7
Rm = 0.2;
B = 0.02;
Rc = 0.02;
Jmax = 1600;
N = 24;
% Current density
J = I/A;
% Voltage losses
Vact = -(z1 + z2*T + z3*T*log(co2) + z4*T*log(I));
Vohm = I*(Rm + Rc);
Vcon = -B*log(1 - J/Jmax);
% Fuel cell voltage
EN=3;
Vfc = N*(EN - Vact - Vohm - Vcon);
% Fuel cell power
P = Vfc .* I;
% Plot
yyaxis left
plot(I,Vfc,'LineWidth',1.5);
ylabel('Fuel Cell Voltage (V)');
xlabel('Fuel Cell Current (A)');
yyaxis right
plot(I,P,'LineWidth',1.5);
ylabel('Fuel Cell Power (W)');
```
 
### Formulae
 
Current density:
 
$$J = \frac{I}{A}$$
 
Activation loss:
 
$$V_{act} = -\left[z_1 + z_2\,T + z_3\,T\ln c_{O_2} + z_4\,T\ln I\right]$$
 
Ohmic loss:
 
$$V_{ohm} = I\,(R_m + R_c)$$
 
Concentration loss:
 
$$V_{con} = -B\ln\left(1 - \frac{J}{J_{max}}\right)$$
 
Fuel cell voltage:
 
$$V_{fc} = N\left(E_N - V_{act} - V_{ohm} - V_{con}\right)$$
 
Fuel cell power:
 
$$P = V_{fc}\,I$$
 
Symbols: $I$ stack current, $A$ cell area, $T$ temperature in kelvin, $z_1$ to $z_4$ activation coefficients, $c_{O_2}$ oxygen concentration (`co2`), $R_m$ membrane resistance, $R_c$ contact resistance, $B$ concentration-loss constant, $J_{max}$ limiting current density, $N$ number of cells, $E_N$ ideal cell voltage (`EN`).
 
---
 
## 5. PV IV Practically Easy (Diode Model)
 
### Code
 
```matlab
clc;
clear;
Rp=400;
Voc=30;
a=2;
I0=5e-05;
Ipv=8;
Vt=1.2;
V=0:0.1:Voc;
I=Ipv-I0.*(exp(V/(a*Vt))-1)-V/Rp;
I(I<0)=0;
P=V.*I;
figure;
plot(V,I,'LineWidth',2)
xlabel('Voltage (V)')
ylabel('Current (A)')
title('I-V Characteristics')
grid on
figure;
plot(V,P,'LineWidth',2)
xlabel('Voltage (V)')
ylabel('Power (W)')
title('P-V Characteristics')
grid on
```
 
### Formulae
 
Current with clamping at zero:
 
$$I = \max\left\{0,\; I_{pv} - I_0\left[\exp\left(\frac{V}{a\,V_t}\right) - 1\right] - \frac{V}{R_p}\right\}$$
 
Power:
 
$$P = V\,I$$
 
Symbols: $I_{pv}$ photocurrent, $I_0$ diode saturation current, $a$ diode ideality factor, $V_t$ thermal voltage, $R_p$ parallel resistance.
 
