

% Program:
% MATLAB INPUT FILE:

clc;
x = input('Enter the number of elements: ');
a = 0;
y = 0;

for e = 1:x
    i = input('Enter the starting node: ');
    j = input('Enter the ending node: ');
    y(e,e) = input('Enter the admittance value: ');
    y(e,e) = 1/y(e,e);
    a(e,i) = 1;
    a(e,j) = -1;
end

ybus = a' * y * a;
ybus

b = input('Enter the number of buses: ');
for i = 2:b
    pg(i) = input('Enter the generator power in MW: ');
    qg(i) = input('Enter the generator power in MVAR: ');
    pL(i) = input('Enter the load in MW: ');
    qL(i) = input('Enter the load in MVAR: ');
    p(i) = pg(i) - pL(i);
    q(i) = qg(i) - qL(i);
end

v(1) = input('Enter the slack bus voltage: ');
acc = input('Enter the acceleration factor: ');

for i = 2:b
    v(i) = input('Enter the voltage: ');
    c(i) = input('Enter the phase angle: ');
    v(i) = v(i) * (exp(1i * (c(i) * (pi / 180))));
    vc(i) = v(i); % Stores initial voltage value
end

n = input('Enter the number of iterations: ');

for p_iter = 1:n
    for i = 2:b
        L(j) = 0;
        for k = 1:b
            if (j ~= k)
                L(j) = L(j) + (ybus(j,k) * v(k));
            end
        end
        s(j) = v(j);
        v(j) = (((p(j) - 1i * q(j)) / conj(v(j))) - L(j)) / ybus(j,j);
        v(j) = s(j) + acc * (v(j) - s(j));
    end
end

for i = 1:b
    disp('bus');
    disp(i);
    disp('voltage: ');
    disp(v(i));
end