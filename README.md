# Power System Programs

## 1. Y-Bus and Z-Bus Formation

```matlab
clc;
clear;

n = input('Enter the number of buses: ');

for i = 1:n
    for j = 1:n
        if(i == j)
            y(i,j) = 0;
            a(i,j) = 0;
        else
            fprintf('Between Bus %d and Bus %d\n', i, j);

            z(i,j) = input('Enter self impedance: ');

            if(z(i,j) == 0)
                y(i,j) = 0;
            else
                y(i,j) = 1/z(i,j);
            end

            a(i,j) = input('Enter half line charging admittance: ');
        end
    end
end

for i = 1:n
    for j = 1:n
        if(i == j)
            R(i,j) = 0;

            for k = 1:n
                R(i,j) = R(i,j) + y(i,k) + a(i,k);
            end

        else
            R(i,j) = -y(i,j);
        end
    end
end

disp('The resultant Y-Bus Matrix is');

R

disp('The resultant Z-Bus Matrix is');

Z = inv(R)
```

## 2. Y-Bus Formation and Gauss-Seidel Load Flow

```matlab
clc;
clear;

e = input('Enter the number of transmission lines: ');
A = [];
Y = zeros(e);
maxBus = 0;

for k = 1:e
    fprintf('\nLine %d\n',k);
    i = input('Enter the starting bus: ');
    j = input('Enter the ending bus: ');
    z = input('Enter the line impedance (R+jX): ');
    maxBus = max([maxBus i j]);
    A(k,maxBus) = 0;
    A(k,i) = 1;
    A(k,j) = -1;
    Y(k,k) = 1/z;
end

Ybus = A' * Y * A;
disp(' ');
disp('Y-Bus Matrix');
disp(Ybus);

b = input('Enter number of buses: ');
baseMVA = input('Enter Base MVA : ');
P = zeros(b,1);
Q = zeros(b,1);

for i = 2:b
    fprintf('\nBus %d\n',i);
    Pg = input('Generator Power (MW): ');
    Qg = input('Generator Reactive Power (MVAR): ');
    Pl = input('Load Power (MW): ');
    Ql = input('Load Reactive Power (MVAR): ');
    P(i) = (Pg-Pl)/baseMVA;
    Q(i) = (Qg-Ql)/baseMVA;
end

V = ones(b,1);
Vm = input('Slack Bus Voltage Magnitude (pu): ');
Va = input('Slack Bus Angle (deg): ');
V(1) = Vm*exp(1i*Va*pi/180);
acc = input('Acceleration Factor (1.0 to 1.8): ');

for i = 2:b
    Vm = input(['Initial Voltage Magnitude of Bus ',num2str(i),' : ']);
    Va = input(['Initial Voltage Angle of Bus ',num2str(i),' : ']);
    V(i) = Vm*exp(1i*Va*pi/180);
end

tol = input('Enter tolerance (e.g. 1e-6): ');
maxIter = input('Enter maximum iterations: ');
fprintf('\nIteration\tMaximum Voltage Change\n');

for iter = 1:maxIter
    Vold = V;
    for i = 2:b
        sumYV = 0;
        for j = 1:b
            if j ~= i
                sumYV = sumYV + Ybus(i,j)*V(j);
            end
        end
        Vnew = (1/Ybus(i,i))*...
               (((P(i)-1i*Q(i))/conj(V(i))) - sumYV);
        V(i) = V(i) + acc*(Vnew - V(i));
    end
    err = max(abs(V-Vold));
    fprintf('%d\t\t%.8f\n',iter,err);
    if err < tol
        fprintf('\nConverged in %d iterations.\n',iter);
        break;
    end
end

disp(' ');
disp('Final Bus Voltages');
for i = 1:b
    fprintf('\nBus %d\n',i);
    fprintf('Voltage Magnitude = %.4f pu\n',abs(V(i)));
    fprintf('Voltage Angle = %.4f degrees\n',angle(V(i))*180/pi);
end
```

## 3. Newton-Raphson Load Flow

```matlab
clc;
clear;

n = input('Enter the number of buses: ');

yb = zeros(n,n);
disp('Enter Y-Bus Matrix (Complex Values)');
for i = 1:n
    for j = i:n
        yb(i,j) = input(['Ybus(',num2str(i),',',num2str(j),') = ']);
        yb(j,i) = yb(i,j);
    end
end

V = zeros(n,1);
delta = zeros(n,1);
Psp = zeros(n,1);
Qsp = zeros(n,1);

for i = 1:n
    V(i) = input(['Voltage Magnitude of Bus ',num2str(i),' = ']);
    delta(i) = deg2rad(input(['Voltage Angle (deg) of Bus ',num2str(i),' = ']));
    Psp(i) = input(['Specified Real Power P(',num2str(i),') = ']);
    Qsp(i) = input(['Specified Reactive Power Q(',num2str(i),') = ']);
end

G = real(yb);
B = imag(yb);
tol = 1e-6;
maxIter = 20;

for iter = 1:maxIter
    P = zeros(n,1);
    Q = zeros(n,1);
    for i = 1:n
        for j = 1:n
            P(i) = P(i) + V(i)*V(j)*( ...
                G(i,j)*cos(delta(i)-delta(j)) + ...
                B(i,j)*sin(delta(i)-delta(j)));
            Q(i) = Q(i) + V(i)*V(j)*( ...
                G(i,j)*sin(delta(i)-delta(j)) - ...
                B(i,j)*cos(delta(i)-delta(j)));
        end
    end

    dP = Psp(2:n) - P(2:n);
    dQ = Qsp(2:n) - Q(2:n);
    mismatch = [dP; dQ];

    fprintf('\nIteration %d\n',iter);
    disp('Mismatch Vector');
    disp(mismatch);

    if max(abs(mismatch)) < tol
        fprintf('\nConverged in %d iterations.\n',iter);
        break;
    end

    J1 = zeros(n-1,n-1);
    J2 = zeros(n-1,n-1);
    J3 = zeros(n-1,n-1);
    J4 = zeros(n-1,n-1);

    for i = 2:n
        for j = 2:n
            if i == j
                J1(i-1,j-1) = -Q(i) - B(i,i)*V(i)^2;
                J2(i-1,j-1) = P(i)/V(i) + G(i,i)*V(i);
                J3(i-1,j-1) = P(i) - G(i,i)*V(i)^2;
                J4(i-1,j-1) = Q(i)/V(i) - B(i,i)*V(i);
            else
                angle = delta(i)-delta(j);
                J1(i-1,j-1) = V(i)*V(j)*( ...
                    G(i,j)*sin(angle) - ...
                    B(i,j)*cos(angle));
                J2(i-1,j-1) = V(i)*( ...
                    G(i,j)*cos(angle) + ...
                    B(i,j)*sin(angle));
                J3(i-1,j-1) = -V(i)*V(j)*( ...
                    G(i,j)*cos(angle) + ...
                    B(i,j)*sin(angle));
                J4(i-1,j-1) = V(i)*( ...
                    G(i,j)*sin(angle) - ...
                    B(i,j)*cos(angle));
            end
        end
    end

    J = [J1 J2;
         J3 J4];

    correction = J \ mismatch;
    dDelta = correction(1:n-1);
    dV = correction(n:end);

    delta(2:n) = delta(2:n) + dDelta;
    V(2:n) = V(2:n) + dV;
end

fprintf('\n----------------------------------\n');
fprintf('FINAL RESULTS\n');
fprintf('----------------------------------\n');
fprintf('\nBus\tVoltage(pu)\tAngle(deg)\n');
for i = 1:n
    fprintf('%d\t%8.4f\t%8.4f\n',i,V(i),rad2deg(delta(i)));
end
```

## 4. Symmetrical and Unsymmetrical Fault Analysis

```matlab
clc;
clear;
close all;

j = 1i;
vpf = 1 + 0j;

Z1 = input('Enter the positive sequence impedance (p.u.) = ');
Z1 = j * Z1;
Z2 = input('Enter the negative sequence impedance (p.u.) = ');
Z2 = j * Z2;
Z0 = input('Enter the zero sequence impedance (p.u.) = ');
Z0 = j * Z0;

Zf = input('Enter the fault impedance (p.u.) = ');
Zf = j * Zf;

Ib = input('Enter the base current (A) = ');

ft = menu('Fault Analysis', ...
    'Three Phase Fault', ...
    'LG Fault', ...
    'LL Fault', ...
    'LLG Fault');

switch ft
    case 1
        If = vpf / (Z1 + Zf);
    case 2
        Ia1 = vpf / (Z1 + Z2 + Z0 + 3*Zf);
        If = 3 * Ia1;
    case 3
        Ia1 = vpf / (Z1 + Z2 + Zf);
        If = -j * sqrt(3) * Ia1;
    case 4
        Ia1 = vpf / (Z1 + (Z2 * (Z0 + 3*Zf)) / (Z2 + Z0 + 3*Zf));
        Ia0 = -Ia1 * (Z2 / (Z2 + Z0 + 3*Zf));
        If = 3 * Ia0;
    otherwise
        error('Invalid choice.');
end

fprintf('\n========== Fault Analysis Result ==========\n');
fprintf('Fault Current (p.u.) = %.4f\n', abs(If));
fprintf('Fault Current (Actual) = %.4f A\n', Ib * abs(If));
fprintf('Fault Current Angle = %.2f degrees\n', rad2deg(angle(If)));
```

## 5. Daily Load Curve and Load Duration Curve (MW)

```matlab
clc;
clear;
close all;

time = [0 6 10 12 16 20 24];
load = [20 25 30 25 35 20];

figure;
stairs(time, [load load(end)], 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Load (MW)');
title('Daily Load Curve');
grid on;

load_sorted = [35 30 25 20];
duration = [4 2 8 10];

cumulative_duration = [0 cumsum(duration)];

figure;
stairs(cumulative_duration, [load_sorted load_sorted(end)], ...
    'LineWidth', 2);
xlabel('Duration (hours)');
ylabel('Load (MW)');
title('Load Duration Curve');
grid on;
```

## 6. Diversity Factor, Load Factor, Load and Load Duration Curves (kW)

```matlab
clc;
clear;
close all;

time = [0 6 8 10 18 24];
load = [100 250 450 300 100];

A = 200;
B = 100;
C = 50;
D = 100;

sum_maximum_demand = A + B + C + D;
station_maximum_demand = max(load);

diversity_factor = sum_maximum_demand / station_maximum_demand;

duration = [6 2 2 8 6];

energy = sum(load .* duration);

average_load = energy / 24;
load_factor = (average_load / station_maximum_demand) * 100;

fprintf('Diversity Factor = %.2f\n', diversity_factor);
fprintf('Units Generated per Day = %.2f kWh/day\n', energy);
fprintf('Load Factor = %.2f %%\n', load_factor);

figure;
stairs(time, [load load(end)], 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Load (kW)');
title('Daily Load Curve');
grid on;

load_sorted = [450 300 250 100];
duration_sorted = [2 8 2 12];

cumulative_duration = [0 cumsum(duration_sorted)];

figure;
stairs(cumulative_duration, [load_sorted load_sorted(end)], ...
    'LineWidth', 2);
xlabel('Duration (hours)');
ylabel('Load (kW)');
title('Load Duration Curve');
grid on;
```

## 7. Diversity Factor and Annual Load Factor

```matlab
clc;
clear;
close all;

Industrial = 1500;
Commercial = 750;
Domestic_power = 100;
Domestic_light = 450;

station_maximum_demand = 2500;
annual_energy = 45e6;

sum_maximum_demand = Industrial + Commercial + ...
                     Domestic_power + Domestic_light;

diversity_factor = sum_maximum_demand / station_maximum_demand;

annual_load_factor = (annual_energy / ...
                     (station_maximum_demand * 8760)) * 100;

fprintf('Diversity Factor = %.2f\n', diversity_factor);
fprintf('Annual Load Factor = %.2f %%\n', annual_load_factor);
```

## 8. Economic Dispatch with Generator Limits

```matlab
clc;
clear;

pd = 925;

a = [0.0045 0.0056 0.0079];
b = [5.2 4.5 5.8];
c = [500 640 820];

B = 0;
A = 0;

pgmax = [450 350 225];
pgmin = [250 200 125];

for i=1:3
    B = B + (b(i)/(2*a(i)));
    A = A + (1/(2*a(i)));
end

lamda = (pd + B) / A;

for i=1:3
    pg(i) = (lamda - b(i)) / (2 * a(i));
end

pg

for i=1:3
    if pg(i)<pgmin(i)
        pgn(i) = pgmin(i);
        k=i;
    elseif pg(i) > pgmax(i)
        pgn(i) = pgmax(i);
        k=i;
    end
end

pgn

pdnw = pd - pgn(k);

Bn = 0;
An = 0;

for i=1:3
    if i~= k
        Bn = Bn + (b(i)/(2*a(i)));
        An = An + (1/(2*a(i)));
    end
end

lamdan = (pdnw + Bn) / An;

for i = 1:3
    if i ~= k
        pgn(i) = (lamdan - b(i)) / (2 * a(i));
    end
end

pgn
```

## 9. Economic Dispatch with Transmission Losses (Loss Coefficients)

```matlab
clc
clear

lambda=18;
n=2;
pd=100;
a=[0.04 0.04];
b=[16 12];
c=[0 0];

B=[0.001 -0.0005
   -0.0005 0.0024];

delp=100;
pg=[25 75];

while delp>0.001

    pgprev=[100 100];
    ppr=[100 100];

    while ppr(1)>1

        for i=1:n

            if i==1
                pg(i)=(lambda-b(i)-(2*lambda*B(i,i+1)*pg(i+1)))/(2*(a(i)+(lambda*B(i,i))));
            else
                pg(i)=(lambda-b(i)-(2*lambda*B(n,n-1)*pg(n-1)))/(2*(a(i)+(lambda*B(i,i))));
            end

            ppr(i)=pgprev(i)-pg(i);
            pgprev(i)=pg(i);

        end

    end

    pg

    PG=0;

    for i=1:n
        PG=PG+pg(i);
    end

    for i=1:n
        loss(i)=(B(i,i)*(pg(i)^2));
    end

    totalloss=2*pg(1)*pg(2)*B(1,2);

    for i=1:n
        totalloss=totalloss+loss(i);
    end

    delp=(pd+totalloss)-PG;

    for i=1:n

        if i==1
            dop(i)=(a(i)+(B(i,i)*b(i))-(2*a(i)*B(i,i+1)*pg(i+1)))/(2*(a(i)+(lambda*B(i,i)))^2);
        else
            dop(i)=(a(i)+(B(i,i)*b(i))-(2*lambda*B(n,n-1)*pg(n-1)))/(2*(a(i)+(lambda*B(i,i)))^2);
        end

    end

    deldeno=0;

    for i=1:n
        deldeno=deldeno+dop(i);
    end

    dellambda=delp/deldeno;
    lambda=lambda+dellambda;

end

pg
PG
totalloss
lambda
```

## 10. Full Load Average Production Cost (Unit Commitment Priority Order)

```matlab
n=3;

a=[0.006 0.01 0.008];

b=[7 8 6];

c=[600 400 500];

pgmax=[400 300 500];

pgmin=[100 50 150];

k=[1.1 1.2 1.0];

for i=1:n

    flapc(i)=(k(i)*(a(i)*pgmax(i)^2+b(i)*pgmax(i)+c(i)))/pgmax(i);

end

flapc

[FLAPC,IX]=sort(flapc);

unit=IX
```
