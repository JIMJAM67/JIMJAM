pd = 925;

a = [0.0045 0.0056 0.0079];
b = [5.2 4.5 5.8];
c = [500 640 820];

B = 0;
A = 0;

pgmax = [450 350 225];
pgmin = [250 200 125];

for i = 1:3
    B = B + (b(i)/(2*a(i)));
    A = A + (1/(2*a(i)));
end

lam = (pd + B)/A;

for i = 1:3
    pg(i) = (lam - b(i))/(2*a(i));
end

disp('Initial Generation:');
disp(pg);

for i = 1:3
    if pg(i) < pgmin(i)
        pgn(i) = pgmin(i);
        k = i;
    elseif pg(i) > pgmax(i)
        pgn(i) = pgmax(i);
        k = i;
    else
        pgn(i) = pg(i);
    end
end

disp('Generation after limits:');
disp(pgn);

pdnw = pd - pgn(k);

Bn = 0;
An = 0;

for i = 1:3
    if i ~= k
        Bn = Bn + (b(i)/(2*a(i)));
        An = An + (1/(2*a(i)));
    end
end

lamnw = (pdnw + Bn)/An;

for i = 1:3
    if i ~= k
        pgn(i) = (lamnw - b(i))/(2*a(i));
    end
end

disp('Final Generation:');
disp(pgn);