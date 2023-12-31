function[Pcyl]=P_isotrans_cyl(C0)

P1111=1/8*(5*C0(1,1)-3*C0(1,2))/C0(1,1)/(C0(1,1)-C0(1,2));
P1122=-1/8*(C0(1,1)+C0(1,2))/C0(1,1)/(C0(1,1)-C0(1,2));
P2222=1/8*(5*C0(1,1)-3*C0(1,2))/C0(1,1)/(C0(1,1)-C0(1,2));
P2323=1/8*1/(0.5*C0(4,4));    % 0.5*C0(4,4)=C0_2323
P1313=1/8*1/(0.5*C0(4,4));    % 0.5*CO(4,4)=C0_2323
P1212=1/8*(3*C0(1,1)-C0(1,2))/C0(1,1)/(C0(1,1)-C0(1,2));

Pcyl=zeros(6,6);

Pcyl(1,1)=P1111;
Pcyl(2,2)=P2222;
Pcyl(1,2)=P1122;
Pcyl(2,1)=Pcyl(1,2);
Pcyl(4,4)=2*P2323;
Pcyl(5,5)=2*P1313;
Pcyl(6,6)=2*P1212;

