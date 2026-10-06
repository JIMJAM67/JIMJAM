clc;
clear;

lambda = 18;
n = 2;
pd = 100;

a = [0.04 0.04];
b = [16 12];

B = [0.001  -0.0005;
    -0.0005  0.0024];

delp = 100;
pg = [25 75];

while abs(delp) > 0.001

    pgprev = [100 100];
    ppr = [100 100];

    while max(abs(ppr)) > 1

        for i = 1:n

            if i == 1
                pg(i) = (lambda - b(i) ...
                    - 2*lambda*B(i,i+1)*pg(i+1)) ...
                    /(2*(a(i) + lambda*B(i,i)));

            else
                pg(i) = (lambda - b(i) ...
                    - 2*lambda*B(i,i-1)*pg(i-1)) ...
                    /(2*(a(i) + lambda*B(i,i)));
            end

            ppr(i) = pgprev(i) - pg(i);
            pgprev(i) = pg(i);

        end

    end

    PG = sum(pg);

    loss = zeros(1,n);
    for i = 1:n
        loss(i) = B(i,i)*pg(i)^2;
    end

    totalloss = 2*pg(1)*pg(2)*B(1,2) + sum(loss);

    delp = (pd + totalloss) - PG;

    for i = 1:n

        if i == 1
            dop(i) = (a(i) + B(i,i)*b(i) ...
                - 2*a(i)*B(i,i+1)*pg(i+1)) ...
                /(2*(a(i) + lambda*B(i,i))^2);

        else
            dop(i) = (a(i) + B(i,i)*b(i) ...
                - 2*lambda*B(i,i-1)*pg(i-1)) ...
                /(2*(a(i) + lambda*B(i,i))^2);
        end

    end

    deldeno = sum(dop);
    dellambda = delp/deldeno;
    lambda = lambda + dellambda;

end

fprintf('\nFinal Results\n');
fprintf('-------------\n');
fprintf('P1 = %.4f MW\n', pg(1));
fprintf('P2 = %.4f MW\n', pg(2));
fprintf('Total Generation (PG) = %.4f MW\n', PG);
fprintf('Transmission Loss = %.4f MW\n', totalloss);
fprintf('Lambda = %.4f\n', lambda);