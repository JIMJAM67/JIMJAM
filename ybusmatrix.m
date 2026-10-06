linedata=[1 6 0.123 0.518 0.01 1;
          1 4 0.08 0.37 0.00 1;
          4 6 0.087 0.407 0.0 1;
          5 2 0.282 0.64 0.0 1;
          2 3 0.723 1.05 0.00 1;
          6 5 0.00 0.3 0.00 1;
          4 3 0.00 0.133 0.00 1];

shuntdata=[1 0.0;
           2 0.0;
           3 0.0;
           4 2.0;
           5 0.0;
           6 2.5];

BaseMVA=100;

n1=linedata(:,1);
nr=linedata(:,2);
R=linedata(:,3);
X=linedata(:,4);
Bc=j*linedata(:,5);
a=linedata(:,6);

busid=shuntdata(:,1);
MVAR=shuntdata(:,2);
sc=j*(MVAR/BaseMVA);

nbr=length(linedata(:,1));
nbus=max(max(n1),max(nr));

Z=R+j*X;
y=1./Z;

% Initialize Ybus
Ybus=zeros(nbus,nbus);

% Off-diagonal elements
for k=1:nbr
    Ybus(n1(k),nr(k)) = Ybus(n1(k),nr(k)) - y(k)/a(k);
    Ybus(nr(k),n1(k)) = Ybus(n1(k),nr(k));
end

% Diagonal elements
for n=1:nbus
    for k=1:nbr
        if n1(k)==n
            Ybus(n,n)=Ybus(n,n)+y(k)/(a(k)^2)+Bc(k);
        end
        if nr(k)==n
            Ybus(n,n)=Ybus(n,n)+y(k)+Bc(k);
        end
    end
end

% Shunt admittances
for n=1:nbus
    idx=find(busid==n);
    if ~isempty(idx)
        Ybus(n,n)=Ybus(n,n)+sc(idx);
    end
end

% Display Ybus and Zbus
disp('Ybus ='); disp(Ybus);
Zbus=inv(Ybus);
disp('Zbus ='); disp(Zbus);
