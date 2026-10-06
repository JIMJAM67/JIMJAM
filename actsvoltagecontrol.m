clc;
clear;

%% SYSTEM DATA
nbus = 3;
nbb = nbus;

linedata = [
    1 2 0.02   0.04   0.00 1;
    1 3 0.01   0.03   0.00 1;
    2 3 0.0125 0.025  0.00 1;
];

busdata = [
    1 1 1.05 0 0.0 0 0 0 0 0;
    2 3 1.00 0 0.0 0 4 2.5 0 0;
    3 2 1.04 0 2.0 0 0 0 0 0;
];

%% INITIALIZATION
bustype = busdata(:,2);
VM = busdata(:,3);
VA = busdata(:,4);

genbus = find((bustype == 2) | (bustype == 1));
loadbus = find(bustype == 3);

itmax = 100;
tol = 1e-12;
nmax = 2*nbus;
basemva = 100;

PGEN  = busdata(:,5)/basemva;
QGEN  = busdata(:,6)/basemva;
PLOAD = busdata(:,7)/basemva;
QLOAD = busdata(:,8)/basemva;

%% SVC DATA
NSVC = 1;
SVCsend = [2];

B = zeros(1,NSVC);
B(1) = 0.01;

BLo = [-0.25];
BHi = [ 0.25];

TarVol = [1.02];
VSta   = [1];

%% Y-BUS FORMATION
j = 1i;

fb = linedata(:,1);
tb = linedata(:,2);

r = linedata(:,3);
x = linedata(:,4);
b = linedata(:,5);
a = linedata(:,6);

b = j*b;

z = r + j*x;
y = 1./z;

nbranch = length(fb);

ybus = zeros(nbus);

% Off-diagonal elements
for k = 1:nbranch
    ybus(fb(k),tb(k)) = ybus(fb(k),tb(k)) - y(k)/a(k);
    ybus(tb(k),fb(k)) = ybus(fb(k),tb(k));
end

% Diagonal elements
for m = 1:nbus
    for n = 1:nbranch

        if fb(n) == m
            ybus(m,m) = ybus(m,m) + y(n)/(a(n)^2) + b(n);

        elseif tb(n) == m
            ybus(m,m) = ybus(m,m) + y(n) + b(n);

        end

    end
end

YR = real(ybus);
YI = imag(ybus);

%% NET INJECTIONS
PNET = (PGEN - PLOAD)';
QNET = (QGEN - QLOAD)';

flag = 0;
it = 1;

%% NEWTON-RAPHSON
while (it < itmax) && (flag == 0)

    %% POWER CALCULATION
    PCAL = zeros(1,nbb);
    QCAL = zeros(1,nbb);

    for ii = 1:nbb

        PSUM = 0;
        QSUM = 0;

        for jj = 1:nbb

            PSUM = PSUM + VM(ii)*VM(jj)* ...
                (YR(ii,jj)*cos(VA(ii)-VA(jj)) + ...
                 YI(ii,jj)*sin(VA(ii)-VA(jj)));

            QSUM = QSUM + VM(ii)*VM(jj)* ...
                (YR(ii,jj)*sin(VA(ii)-VA(jj)) - ...
                 YI(ii,jj)*cos(VA(ii)-VA(jj)));

        end

        PCAL(ii) = PSUM;
        QCAL(ii) = QSUM;
    end

    %% SVC REACTIVE INJECTION
    for ii = 1:NSVC
        QCAL(SVCsend(ii)) = ...
            QCAL(SVCsend(ii)) - VM(SVCsend(ii))^2*B(ii);
    end

    %% POWER MISMATCHES
    DP = PNET - PCAL;
    DQ = QNET - QCAL;

    for ii = 1:nbb

        if bustype(ii) == 1
            DP(ii) = 0;
            DQ(ii) = 0;

        elseif bustype(ii) == 2
            DQ(ii) = 0;
        end

    end

    DPQ = zeros(nmax,1);

    kk = 1;
    for ii = 1:nbb

        DPQ(kk)   = DP(ii);
        DPQ(kk+1) = DQ(ii);

        kk = kk + 2;

    end

    %% CONVERGENCE CHECK
    if max(abs(DPQ)) < tol
        flag = 1;
        break;
    end

    %% JACOBIAN
    JAC = zeros(nmax,nmax);

    iii = 1;

    for ii = 1:nbb

        jjj = 1;

        for jj = 1:nbb

            if ii == jj

                JAC(iii,jjj) = -QCAL(ii) - VM(ii)^2*YI(ii,ii);
                JAC(iii,jjj+1) = PCAL(ii) + VM(ii)^2*YR(ii,ii);

                JAC(iii+1,jjj) = PCAL(ii) - VM(ii)^2*YR(ii,ii);
                JAC(iii+1,jjj+1) = QCAL(ii) - VM(ii)^2*YI(ii,ii);

            else

                JAC(iii,jjj) = VM(ii)*VM(jj)* ...
                    (YR(ii,jj)*sin(VA(ii)-VA(jj)) - ...
                     YI(ii,jj)*cos(VA(ii)-VA(jj)));

                JAC(iii+1,jjj) = -VM(ii)*VM(jj)* ...
                    (YI(ii,jj)*sin(VA(ii)-VA(jj)) + ...
                     YR(ii,jj)*cos(VA(ii)-VA(jj)));

                JAC(iii,jjj+1) = -JAC(iii+1,jjj);
                JAC(iii+1,jjj+1) = JAC(iii,jjj);

            end

            jjj = jjj + 2;

        end

        iii = iii + 2;

    end

    %% SLACK & PV MODIFICATION

    for kk = 1:nbb

        if bustype(kk) == 1

            row = kk*2 - 1;

            for jj = 1:2*nbb

                if row == jj
                    JAC(row,row) = 1;
                else
                    JAC(row,jj) = 0;
                    JAC(jj,row) = 0;
                end

            end

        end

        if (bustype(kk) == 1) || (bustype(kk) == 2)

            row = kk*2;

            for jj = 1:2*nbb

                if row == jj
                    JAC(row,row) = 1;
                else
                    JAC(row,jj) = 0;
                    JAC(jj,row) = 0;
                end

            end

        end

    end

    %% SVC JACOBIAN UPDATE
    for ii = 1:NSVC

        if VSta(ii) == 1

            JAC(:,2*SVCsend(ii)) = 0;

            JAC(2*SVCsend(ii)-1,2*SVCsend(ii)-1) = ...
                JAC(2*SVCsend(ii)-1,2*SVCsend(ii)-1) ...
                - VM(SVCsend(ii))^2*B(ii);

            JAC(2*SVCsend(ii),2*SVCsend(ii)) = ...
                -VM(SVCsend(ii))^2*B(ii);
        end

    end

    %% SOLVE
    D = JAC \ DPQ;

    %% UPDATE STATES
    iii = 1;

    for ii = 1:nbb

        VA(ii) = VA(ii) + D(iii);
        VM(ii) = VM(ii) + D(iii+1)*VM(ii);

        iii = iii + 2;

    end

    %% UPDATE SVC
    for ii = 1:NSVC

        if VSta(ii) == 1

            VM(SVCsend(ii)) = TarVol(ii);

            value = B(ii)*D(2*SVCsend(ii));
            value2 = D(2*SVCsend(ii));

            if value > 0.1
                value2 = 0.1/B(ii);
            elseif value < -0.1
                value2 = -0.1/B(ii);
            end

            B(ii) = B(ii) + B(ii)*value2;

        end

    end

    for ii = 1:NSVC

        if B(ii) > BHi(ii)
            B(ii) = BHi(ii);
        elseif B(ii) < BLo(ii)
            B(ii) = BLo(ii);
        end

    end

    it = it + 1;

end

%% RESULTS
VA1 = VA*180/pi;

disp('------------------------------');
disp('|  Bus  |    V    |  Angle   |');
disp('|  No   |   pu    | Degree   |');
disp('------------------------------');

for m = 1:nbb
    fprintf('%4d %12.6f %12.6f\n',m,VM(m),VA1(m));
end

disp('------------------------------');
disp(['Iterations = ',num2str(it)]);
disp(['SVC B = ',num2str(B)]);
