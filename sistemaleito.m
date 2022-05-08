function F = sistemaleito(x,Qe,QF,N1F,N2F,N3F,N4F,N5F,N6F,Te,uo,umf,Ab,Z,V,teta1,teta2,teta3,teta4,teta5,teta6,ef,Kr,Ka,Cc,P,R,delta)
    F = [
        x(1)-(Qe/QF)*N1F-(uo-umf)*Ab*(N1F/QF-x(1)/Qe)*(1-exp(-teta1*Z))+V*(1-delta)*(1-ef)*((Kr*Ka*Cc*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)))/((R*Te)+(Ka*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)));
        x(2)-(Qe/QF)*N2F-(uo-umf)*Ab*(N2F/QF-x(2)/Qe)*(1-exp(-teta2*Z))+2* V*(1-delta)*(1-ef)*((Kr*Ka*Cc*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)))/((R*Te)+(Ka*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)));
        x(3)-(Qe/QF)*N3F-(uo-umf)*Ab*(N3F/QF-x(3)/Qe)*(1-exp(-teta3*Z))-V*(1-delta)*(1-ef)*((Kr*Ka*Cc*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)))/((R*Te)+(Ka*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)));
        x(4)-(Qe/QF)*N4F-(uo-umf)*Ab*(N4F/QF-x(4)/Qe)*(1-exp(-teta4*Z))+0.5*V*(1-delta)*(1-ef)*((Kr*Ka*Cc*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)))/((R*Te)+(Ka*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)));
        x(5)-(Qe/QF)*N5F-(uo-umf)*Ab*(N5F/QF-x(5)/Qe)*(1-exp(-teta5*Z))-V*(1-delta)*(1-ef)*((Kr*Ka*Cc*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)))/((R*Te)+(Ka*x(1)*P)/(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)));
        x(6)-(Qe/QF)*N6F-(uo-umf)*Ab*(N6F/QF-x(6)/Qe)*(1-exp(-teta6*Z))-V*(1-delta)*(1-ef)*0;
    ];
end