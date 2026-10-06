If=0;
j=sqrt(-1);
vpf=1+j*0;
Z1=input('enter the positive sequence impedance of the system in p.u. = ');
Z1=j*Z1;
Z2=input('enter the negative sequence impedance of the system in p.u. = ');
Z2=j*Z2;
Z0=input('enter the zero sequence impedance of the systeminp.u. = ');
Z0=j*Z0;
Zf=input('Enter the fault impedanceinp.u. = ');
Zf=j*Zf;
Ib=input('Enter the base current = ');
ft=menu('option','Three Phase Fault','LGFault','LLFault','LLG Fault')
if(ft==1)
        If=vpf/Z1+Zf;
If_magnitude_in_perunit=abs(If)
If_magnitude_in_acutualvalue=Ib*abs(If)
If_angle=(180/pi)*angle(If)
end
if(ft==2)
        Ia1=vpf/(Z1+Z2+Z0+3*Zf);
       If=3*Ia1;
If_magnitude_in_perunit=abs(If)
If_magnitude_in_acutualvalue=Ib*abs(If)
If_angle=(180/pi)*angle(If)
end
if(ft==3)
        Ia1=vpf/(Z1+Z2+Zf);
        If=-j*sqrt(3)*Ia1;
If_magnitude_in_perunit=abs(If)
If_magnitude_in_acutualvalue=Ib*abs(If)
If_angle=(180/pi)*angle(If)
end
if(ft==4)
        Ia1=vpf/(Z1+((Z2*(Z0+3*Zf))/(Z2+Z0+3*Zf)));
       Ia0=-Ia1*(Z2/(Z2+Z0+3*Zf));
       If=3*Ia0;
If_magnitude_in_perunit=abs(If)
If_magnitude_in_acutualvalue=Ib*abs(If)
If_angle=(180/pi)*angle(If)
end


