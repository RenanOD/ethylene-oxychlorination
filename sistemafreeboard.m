function dN = sistemafreeboard(Z,N,ef,A,a,Kr,Ka,Cc,P,R,T,N6)
    dN = zeros(5,1);
    dN(1)=(1-ef)*A*exp(-a*Z)*(((-1*Kr*Ka*Cc*N(1)*P)/(N(1)+N(2)+N(3)+N(4)+N(5)+N6))/(((R*T)+(Ka*N(1)*P))/(N(1)+N(2)+N(3)+N(4)+N(5)+N6)));
    dN(2)=(1-ef)*A*exp(-a*Z)*(((-2*Kr*Ka*Cc*N(1)*P)/(N(1)+N(2)+N(3)+N(4)+N(5)+N6))/(((R*T)+(Ka*N(1)*P))/(N(1)+N(2)+N(3)+N(4)+N(5)+N6)));
    dN(3)=(1-ef)*A*exp(-a*Z)*(((Kr*Ka*Cc*N(1)*P)/(N(1)+N(2)+N(3)+N(4)+N(5)+N6))/(((R*T)+(Ka*N(1)*P))/(N(1)+N(2)+N(3)+N(4)+N(5)+N6)));
    dN(4)=(1-ef)*A*exp(-a*Z)*(((-0.5*Kr*Ka*Cc*N(1)*P)/(N(1)+N(2)+N(3)+N(4)+N(5)+N6))/(((R*T)+(Ka*N(1)*P))/(N(1)+N(2)+N(3)+N(4)+N(5)+N6)));
    dN(5)=(1-ef)*A*exp(-a*Z)*(((Kr*Ka*Cc*N(1)*P)/(N(1)+N(2)+N(3)+N(4)+N(5)+N6))/(((R*T)+(Ka*N(1)*P))/(N(1)+N(2)+N(3)+N(4)+N(5)+N6)));
end