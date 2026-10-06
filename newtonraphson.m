a% ========================================================================
% Newton-Raphson Load Flow using Y-bus (rectangular form G, B)
% Bus 1 is treated as the slack bus (its angle and magnitude are fixed).
% ========================================================================

clear all;
clc;

n = input('Enter the number of buses: ');

% ---- Read Y-bus matrix ----
for i = 1:n
    for j = 1:n
        Ybus(i,j) = input(['Enter the Y bus matrix element ', num2str(i), num2str(j), ': ']);
    end
end

% ---- Read initial bus data ----
for i = 1:n
    mag(i) = input(['Enter the voltage magnitude of bus ', num2str(i), ': ']);
    th(i)  = input(['Enter the angle of bus ',             num2str(i), ': ']);
    aep(i) = input(['Enter the real power of bus ',        num2str(i), ': ']);
    acq(i) = input(['Enter the reactive power of bus ',    num2str(i), ': ']);
end

my = abs(Ybus);      % Y-bus magnitude (not used directly below, kept for reference)
an = angle(Ybus);    % Y-bus angle     (not used directly below, kept for reference)
g  = real(Ybus);     % conductance matrix G
b  = imag(Ybus);     % susceptance matrix B

tol     = 1e-6;
maxiter = 50;
iter    = 0;

while true
    iter = iter + 1;

    % ---- Compute injected P and Q at every bus ----
    for i = 1:n
        pe(i) = 0;
        qu(i) = 0;
        for j = 1:n
            pe(i) = pe(i) + mag(i)*mag(j)*( g(i,j)*cos(th(i)-th(j)) + b(i,j)*sin(th(i)-th(j)) );
            qu(i) = qu(i) + mag(i)*mag(j)*( g(i,j)*sin(th(i)-th(j)) - b(i,j)*cos(th(i)-th(j)) );
        end
    end

    % ---- Build the Jacobian (sizes n x n, slack row/col removed later) ----
    for i = 2:n
        for j = 2:n
            if i ~= j
                J1(i,j) =  mag(i)*mag(j)*( g(i,j)*sin(th(i)-th(j)) - b(i,j)*cos(th(i)-th(j)) );
                J3(i,j) = -mag(i)*mag(j)*( g(i,j)*cos(th(i)-th(j)) + b(i,j)*sin(th(i)-th(j)) );
                J2(i,j) = -J3(i,j);
                J4(i,j) =  J1(i,j);
            else
                J1(i,j) = -qu(i) - b(i,i)*mag(i)^2;
                J2(i,j) =  pe(i) + g(i,i)*mag(i)^2;
                J3(i,j) =  pe(i) - g(i,i)*mag(i)^2;
                J4(i,j) =  qu(i) - b(i,i)*mag(i)^2;
            end
        end
    end

    % ---- Assemble reduced Jacobian (slack bus excluded) ----
    Ja1 = J1(2:n, 2:n);
    Ja2 = J2(2:n, 2:n);
    Ja3 = J3(2:n, 2:n);
    Ja4 = J4(2:n, 2:n);
    Jacob = [Ja1 Ja2; Ja3 Ja4];

    % ---- Power mismatches ----
    delp = (aep(2:n) - pe(2:n))';
    delq = (acq(2:n) - qu(2:n))';

    mismatch = [delp; delq];

    % ---- Convergence check ----
    if max(abs(mismatch)) < tol || iter >= maxiter
        break;
    end

    % ---- Solve for corrections ----
    corr = Jacob \ mismatch;     % corr(1:n-1) = d(theta), corr(n:2*(n-1)) = dV/V

    chth  = zeros(1,n);
    chmag = zeros(1,n);
    chth(2:n)  = corr(1:n-1)';
    chmag(2:n) = corr(n:2*(n-1))';

    % ---- Undo the ΔV/V normalization to get actual ΔV ----
    for i = 2:n
        chmag(i) = chmag(i) * mag(i);
    end

    % ---- Update state ----
    th  = th  + chth;
    mag = mag + chmag;
end

disp(['Converged in ', num2str(iter), ' iterations']);
disp('The voltage magnitudes are:');
disp(mag);
disp('The phase angles (rad) are:');
disp(th);
disp('The bus injected real powers are:');
disp(pe);
disp('The bus injected reactive powers are:');
disp(qu);