% ARQUIVO PRINCIPAL
%
% Objetivo: Simular reator de leito fluidizado
%
%------------------------------------------------------------------------
% Definição dos componentes
%------------------------------------------------------------------------
% 1 - Etileno
% 2 - HCl
% 3 - H2O
% 4 - O2
% 5 - EDC
% 6 - N2
%------------------------------------------------------------------------
% Definição das Constantes e Propriedades física
%------------------------------------------------------------------------
% NiF - Vazões molares na entrada do reator (kgmol/s)
% Nie - Vazões molares na fase emulsão (kgmol/s)
% Nib - Vazões molares na fase bolha (kgmol/s)
% Nout(i)- Vazões molares na saída (kgmol/s)
% Te  - Temperatura na fase emulsão (K)
% Tref- Temperatura de referência (K)
% TMF - Temperatura média da corrente de alimentação (K)
% Tw  - Temperatura 
% Tc  - Temperaturas críticas dos componentes (K)
% Pc  - Pressões críticas dos componentes (Pa)
% Vc  - Volume crítico dos componentes (m3/kgmol)
% Zc  - Fator de compressibilidade crítico (-)
% R   - Constante universal dos gases 8314 (Pa.m3/kgmol.K)
% Ka  - Constante de equilíbrio da adsorção (m3/kgmol)
% Cc  - Concentração de cloreto de cobre no catalisador (kgmol/m3_cat)
% Psi - Esfericidade da partícula de catalisador (-)
% g   - Aceleração da gravidade (m/s2)
% m   - Vazões mássicas dos reagente (kg/s)
% M   - Pesos moleculares dos componentes (kgmol)
% DENs- Massa específica do catalisador (kg/m3)
% Dp  - Momento dipolar (-)
% At  - Área da seção transversal do leito (m2)
% u0  - Velocidade superficial do gás (m/s)
% P   - Pressão no interior do reator (Pa)
% Aw  - Área da serpentina para troca térmica (m2)
% kg  - Condutividade térmica do fluido na emulsão (J/m.s.K)
% hw  - Coeficiente de película fase emulsão (J/m2.s.K)
% hi  - Coeficiente de película no interior do tubo (J/m2.s.K)
% r1  - Velocidade da reação principal por volume de cat (kgmol/s.m3_cat)
% L   - Comprimento do tubo de troca térmica (m)
% ri  - Raio interno do tubo de troca térmica (m)
% re  - Raio externo do tubo de troca térmica (m)
% Kt  - Condutividade térmica do tubo (J/s.m.K)
% U   - Coeficiente global de troca térmica (J/m2.s.K)
% V	- Volumes atômicos difusionais (Angstron);
% W   - Massa de catalisador (kg)
%------------------------------------------------------------------------
% Definição das variáveis Físico-químicas e Termodinâmicas
%------------------------------------------------------------------------
% yiF   - Fração molar dos componentes na entrada (-)
% yie   - Fração molar dos componentes na emulsão (-)
% DENgF - Densidade média dos gases na entrada (kg/m3)
% DENgo - Densidade média dos gases na saída (kg/m3)
% CPgF  - Capacidade calorífica molar média do gás na entrada (J/kgmol.K)
% CPgo  - Capacidade calorífica molar media do gás na saída (J/kgmol.K)
% CPgl  - Capacidade calirífica molar média do gás no leito (J/kgmol.K)
% mi    - Viscosidade do gás no leito (Pa.s)
% DH    - Calor da reação (J/kgmol.s)
% Di    - Difusividade dos componentes na mistura (m2/s)
%------------------------------------------------------------------------
% Definição dos Parâmetros hidrodinâmicos e Coeficientes de transferência
%------------------------------------------------------------------------
% Qo    - Vazão volumétrica na saída da fase emulsão (m3/s)
% QF    - Vazão volumétrica de alimentação (m3/s)
% Ab    - Área da seção transversal da fase bolha (m2)
% At    - Área da seção transversal do reator (m2)
% teta - Constante da equação do balanço de massa
% beta  - Variável beta
% umf   - Velocidade mínima de fluidização (m/s)
% db    - Diâmetro das bolhas (m)
% ub    - Velocidade de ascensão das bolhas (m/s)
% Kbc   - Coef de transferência de massa bolha/suvem (m/s)
% Kce   - Coef de transferência de massa nuvem/emulsão (m/s)
% Kr    - Constante cinética da reação de oxicloração (1/s)
% Hbe   - Coef de transferência de calor bolha/emulsão (J/m3.s.K)
%
%--------------------------- INÍCIO DO PROGRAMA -------------------------
%------------------------------------------------------------------------
% Estimativas iniciais
%------------------------------------------------------------------------
N1e=1.49*10^-3; N2e=0.0029; N3e=0.0538; N4e=0.0137; N5e=0.0538; N6e=0.032; 
Nout(1)=1.49*10^-3; Nout(2)=0.0029; Nout(3)=0.0538;
Nout(4)=0.0137; Nout(5)=0.0538; Nout(6)=0.032;
Te=543; TMF=543; R=8314; PMF=840000;
m1=0.42951; m2=1.11618; mar=1.35626;
M=[28.06,36.46,18.02,32,98.96,28.01];
N1F=m1/M(1); N2F=m2/M(2); Nar=mar/28.8479; mN2entrada=Nar*0.79*M(6);
Qentrada=(N1F+N2F+Nar)*R*TMF/PMF;
DENentrada=(PMF/(R*TMF))*(N1F*M(1)+N2F*M(2)+Nar*28.8479)/(N1F+N2F+Nar);
%------------------------------------------------------------------------
% Contadores
%------------------------------------------------------------------------
conta_energia=0; conta_massa=0; loop_energia=1;
while loop_energia==1           % While será fechado no final do programa
%------------------------------------------------------------------------
% Constantes e parâmetros de entrada
%------------------------------------------------------------------------
CP=[4.221   -8.782  5.795   -6.729  2.511;
    3.827   -2.936  0.879   -1.031  0.439;
    4.395   -4.186  1.405   -1.564  0.632;
    3.63    -1.794  0.658   -0.601  0.179;
    2.99    23.197  -0.404  -1.133  0.617;
    3.539   -0.261  0.007   0.157   -0.099];
Tref=298;
m1=0.42951; m2=1.11618; mar=1.35626;
M=[28.06,36.46,18.02,32,98.96,28.01];
V=[41.04,23.31,10.73,12.22,83.04,9.08];                
Tc=[282.34,324.61,647.14,154.58,561.00,126.2];
Pc=[5041000,8310000,22064000,5043000,5400000,3394387];
Zc=[0.282,0.249,0.229,0.288,0.255,0.290];
Vc=[0.131,0.081,0.0559,0.0734,0.220,0.0895];
Dp=[0,1.1,1.8,0,1.8,0];                                %Momentos dipolar
k=[0,0,0,0,0,0];                                       %Fator correção cond. térmica (adm)
omega=[0.087,0,0.344,0,0,0.037];                       %Fator acênctrico
R=8314;                                                %Constante dos gases (J/kgmol.K)
Ka=630;                                                %Constante de equilíbrio de adsorção (m3/kgmol)
g=9.81;                                                %Aceleração da gravidade (m/s2)
DENs=1369;                                             %Densidade da partícula de catalisador (kg/m3)
DI=2.5;                                                %Diâmetro interno do reator(m)
P=735000;                                              %Pressão de operação no reator (Pa)
Tw=313;                                                %Temperatura da água de resfriamento (K)
dp=190*10^-6;                                          %Diâmetro da partícula (m) 
EAP=6/(dp*DENs);                                       %Área superficia da partícula (m2/kg)
W=12750;                                                %Massa de catalisador no leito (Kg)
dorificio=1/(50*pi^0.5);                               %Tal que norificios/m2 = 1000
Ao=pi*dorificio^2/4;                                  %Área de cada orifício do distribuidor;
norificios=DI^2/(25*dorificio^2);
pitch=dorificio/(0.04/0.9)^0.5;                  %Pitch triangular
%------------------------------------------------------------------------
% Cálculo das frações molares na entrada
%------------------------------------------------------------------------
N1F=m1/M(1);            N2F=m2/M(2);    N3F=0;
N4F=0.21*mar/28.8479;     N5F=0;          N6F=0.79*mar/28.8479;
Nar=mar/28.8479;
NTF=N1F+N2F+Nar;
y1F=N1F/NTF; y2F=N2F/NTF; y3F=N3F/NTF;
y4F=N4F/NTF; y5F=N5F/NTF; y6F=N6F/NTF;
yar=Nar/NTF;
%------------------------------------------------------------------------
% Cálculo das frações molares na emulsão
%------------------------------------------------------------------------
NTe=N1e+N2e+N3e+N4e+N5e+N6e;
y1e=N1e/NTe; y2e=N2e/NTe; y3e=N3e/NTe; y4e=N4e/NTe; y5e=N5e/NTe; y6e=N6e/NTe;
%------------------------------------------------------------------------
% Cálculo das frações molares na saída
%------------------------------------------------------------------------
NoutTot=Nout(1)+Nout(2)+Nout(3)+Nout(4)+Nout(5)+Nout(6);
y1o=Nout(1)/NoutTot; y2o=Nout(2)/NoutTot; y3o=Nout(3)/NoutTot;
y4o=Nout(4)/NoutTot; y5o=Nout(5)/NoutTot; y6o=Nout(6)/NoutTot;
%------------------------------------------------------------------------
% Cálculo dos parâmetros na entrada
%------------------------------------------------------------------------
TMF=543;
PMF=840000;
QF=(Nar+N1F+N2F)*R*TMF/PMF;
MMF=y1F*M(1)+y2F*M(2)+yar*28.8479;
MMe=y1e*M(1)+y2e*M(2)+y3e*M(3)+y4e*M(4)+y5e*M(5)+y6e*M(6);
MMo=y1e*M(1)+y2e*M(2)+y3e*M(3)+y4e*M(4)+y5e*M(5)+y6e*M(6);
%------------------------------------------------------------------------
% Cálculo do fator de compressibilidade
%------------------------------------------------------------------------
for i=1:6
    a(i)=0.42748*(R^2*Tc(i)^2.5/Pc(i));
    b(i)=0.08664*(R*Tc(i)/Pc(i));
end
bmF=y1F*b(1)+y2F*b(2)+y3F*b(3)+y4F*b(4)+y5F*b(5)+y6F*b(6);
bme=y1e*b(1)+y2e*b(2)+y3e*b(3)+y4e*b(4)+y5e*b(5)+y5e*b(6);
amF=(y1F*(a(1)^0.5)+y2F*(a(2)^0.5)+y3F*(a(3)^0.5)+y4F*(a(4)^0.5)+y5F*(a(5)^0.5)+y6F*(a(6)^0.5))^2;
ame=(y1e*(a(1)^0.5)+y2e*(a(2)^0.5)+y3e*(a(3)^0.5)+y4e*(a(4)^0.5)+y5e*(a(5)^0.5)+y6e*(a(6)^0.5))^2;
VidF=R*TMF/PMF;
Vide=R*Te/P;
VcalcF=@(VF)(PMF+amF/(TMF^0.5*(VF+bmF)*VF))*(VF-bmF)-R*TMF;
Vcalce=@(Ve)(P+ame/(Te^0.5*(Ve+bme)*Ve))*(Ve-bme)-R*Te;
VrF=fzero(VcalcF,VidF);
Vre=fzero(Vcalce,Vide);
ZF=VrF/VidF;
Ze=Vre/Vide;
%------------------------------------------------------------------------
% Cálculo das densidades na entrada, leito e saída
%------------------------------------------------------------------------
DENgF=PMF*MMF/(ZF*R*TMF);
DENgl=P*MMe/(Ze*R*Te);
DENgo=P*MMo/(Ze*R*Te);
%------------------------------------------------------------------------
% Cálculo da difusividade dos componentes na mistura
%------------------------------------------------------------------------
for i=1:6
    for j=1:6
        D(i,j)=(0.001*Te^1.75*(1/M(1,i)+1/M(1,j))^0.5)/(P*9.869*10^-6*((V(1,i)^0.33+V(1,j)^0.33))^2);
    end
end
Di(1)=(1-y1e)*10^-4/((y2e/D(1,2))+(y3e/D(1,3))+(y4e/D(1,4))+(y5e/D(1,5))+(y6e/D(1,6)));
Di(2)=(1-y2e)*10^-4/((y1e/D(2,1))+(y3e/D(2,3))+(y4e/D(2,4))+(y5e/D(2,5))+(y6e/D(2,6)));
Di(3)=(1-y3e)*10^-4/((y1e/D(3,1))+(y2e/D(3,2))+(y4e/D(3,4))+(y5e/D(3,5))+(y6e/D(3,6)));
Di(4)=(1-y4e)*10^-4/((y1e/D(4,1))+(y2e/D(4,2))+(y3e/D(4,3))+(y5e/D(4,5))+(y6e/D(4,6)));
Di(5)=(1-y5e)*10^-4/((y1e/D(5,1))+(y2e/D(5,2))+(y3e/D(5,3))+(y4e/D(5,4))+(y6e/D(5,6)));
Di(6)=(1-y6e)*10^-4/((y1e/D(6,1))+(y2e/D(6,2))+(y3e/D(6,3))+(y4e/D(6,4))+(y5e/D(6,5)));
%------------------------------------------------------------------------
% Cálculo das capacidades caloríficas na entrada, leito e saída
%------------------------------------------------------------------------
for i=1:6
    CPF(i)=(CP(i,1)+CP(i,2)*10^-3*TMF+CP(i,3)*10^-5*TMF^2+CP(i,4)*10^-8*TMF^3+CP(i,5)*10^-11*TMF^4)*R;
end
for i=1:6
    CPo(i)=(CP(i,1)+CP(i,2)*10^-3*512.5+CP(i,3)*10^-5*512.5^2+CP(i,4)*10^-8*512.5^3+CP(i,5)*10^-11*512.5^4)*R;
end
for i=1:6
    CPl(i)=(CP(i,1)+CP(i,2)*10^-3*Te+CP(i,3)*10^-5*Te^2+CP(i,4)*10^-8*Te^3+CP(i,5)*10^-11*Te^4)*R;
end
CPFar=(3.355+0.575*10^-3*TMF-0.016*10^5*TMF^-2)*R;
CPgF=(CPF(1)*y1F+CPF(2)*y2F+CPFar*yar);
CPgo=(CPl(1)*y1e+CPl(2)*y2e+CPl(3)*y3e+CPl(4)*y4e+CPl(5)*y5e+CPl(6)*y6e);
CPgl=(CPo(1)*y1o+CPo(2)*y2o+CPo(3)*y3o+CPo(4)*y4o+CPo(5)*y5o+CPo(6)*y6o);
%------------------------------------------------------------------------
% Cálculo da viscosidade no leito
%------------------------------------------------------------------------
Mm=y1e*M(1)+y2e*M(2)+y3e*M(3)+y4e*M(4)+y5e*M(5)+y6e*M(6);
Tcm=y1e*Tc(1)+y2e*Tc(2)+y3e*Tc(3)+y4e*Tc(4)+y5e*Tc(5)+y6e*Tc(6);
Zcm=y1e*Zc(1)+y2e*Zc(2)+y3e*Zc(3)+y4e*Zc(4)+y5e*Zc(5)+y6e*Zc(6);
Vcm=(y1e*Vc(1)+y2e*Vc(2)+y3e*Vc(3)+y4e*Vc(4)+y5e*Vc(5)+y6e*Vc(6))*1000;
Pcm=R*10^-2*Tcm*(Zcm/Vcm);
Trm=Te/Tcm;
Prm=P*10^-5/Pcm;
zeta=0.176*(Tcm/(Mm^3*Pcm^4))^(1/6);
for i=1:6
    Tr(i)=Te/Tc(i);
end
for i=1:6
    Dpr(i)=(52.46*Dp(i)^2*Pc(i))/Tc(i)^2;
end
for i=1:6
    if Dpr(i)>=0 & Dpr(i)<0.022
        Fpo(i)=1;
    elseif Dpr(i)>=0.022 & Dpr(i)<0.075;
        Fpo(i)=1+30.55*(0.292-Zc(i))^1.72;
    elseif Dpr(i)>=0.075;
        Fpo(i)=1+30.55*(0.292-Zc(i))^1.72*abs((0.96+0.1*(Tr(i)-0.7)));
    end
end
Fpom=y1e*Fpo(1)+y2e*Fpo(2)+y3e*Fpo(3)+y4e*Fpo(4)+y5e*Fpo(5)+y6e*Fpo(6);
Z1=((0.807*Trm^0.618-0.357*exp(-0.449*Trm)+0.340*exp(-4.058*Trm)+0.018))*Fpom;
Visco=Z1/zeta;
mi=Visco*10^-7;
%------------------------------------------------------------------------
% Cálculo da constante de Arrhenius
%------------------------------------------------------------------------
Kr=269*exp(-37.8/(R*Te*10^-6));
%------------------------------------------------------------------------
% Cálculo da condutividade térmica do leito
%------------------------------------------------------------------------
y=[y1e,y2e,y3e,y4e,y5e,y6e];
for i=1:6
    p(i)=0.809*(Vc(i)*10^3)^(1/3);
end
for i=1:6
    for j=1:6
        pp(i,j)=(p(i)*p(j))^0.5;
    end
end
for i=1:6
    for j=1:6
        t(i,j)=y(i)*y(j)*pp(i,j)^3;
    end
end
pmtres=sum(t);
pmtres=sum(pmtres);
pm=(pmtres)^(1/3);
for i=1:6
    sk1(i)=Tc(i)/1.2593;
end
for i=1:6
    for j=1:6
        sk(i,j)=(sk1(i)*sk1(j))^0.5;
    end
end
for i=1:6
    for j=1:6
        t(i,j)=y(i)*y(j)*pp(i,j)^3*sk(i,j);
    end
end
c=sum(t);
c=sum(c);
skm=c/pmtres;
for i=1:6
    for j=1:6
        M1(i,j)=(2*M(i)*M(j))/(M(i)+M(j));
    end
end
for i=1:6
    for j=1:6
        t(i,j)=y(i)*y(j)*pp(i,j)^2*M1(i,j)^0.5*sk(i,j);
    end
end
d=sum(t);
d=sum(d);
e=skm*pm^2;
Mmm=(d/e)^2;
for i=1:6
    for j=1:6
        omega1(i,j)=(omega(i)+omega(j))/2;
    end
end
for i=1:6
    for j=1:6
        omega1(i,j)=(omega(i)+omega(j))/2;
    end
end
for i=1:6
    for j=1:6
        t(i,j)=y(i)*y(j)*omega1(i,j)*pp(i,j)^3;
    end
end
f=sum(t);
f=sum(f);
omegam=f/pm^3;
for i=1:6
    for j=1:6
        t(i,j)=(y(i)*y(j)*Dp(i)^2*Dp(j)^2/pp(i,j)^3);
    end
end
h=sum(t);
h=sum(h);
Dpm=(pm^3*h)^0.25;
Vcmm=(pm/0.809)^3;
Tcmm=1.2593*skm;
Dprm=131.3*Dpm/(Vcmm*Tcmm)^0.5;
for i=1:6
    for j=1:6
        k(i,j)=(k(i)*k(j))^0.5;
    end
end
for i=1:6
    for j=1:6
        t(i,j)=y(i)*y(j)*k(i,j);
    end
end
km=sum(t);
km=sum(km);
Fcm=1-0.275*omegam+0.059035*Dprm^4;
Tam=Te/skm;
ww=1.16145*Tam^-0.14874+0.52487*exp(-0.77320*Tam)+2.16178*exp(-2.43787*Tam);
nm=26.69*Fcm*(Mmm*Te)^0.5/(pm^2*ww);

for i=1:6
    CV(i)=(CPl(i)*10^-3)-R*10^-3;
end
CVm=y1e*CV(1)+y2e*CV(2)+y3e*CV(3)+y4e*CV(4)+y5e*CV(5)+y6e*CV(6);
Trmm=Te/Tcmm;
Mlm=Mmm/1000;
alfam=(CVm/(R*10^-3))-(2/3);
betam=0.7862-0.7109*omegam+1.316*omegam^2;
Zm=2+10.5*Trmm^2;
td=1+alfam*((0.215+0.28288*alfam-1.061*betam+0.26665*Zm)/(0.6366+betam*Zm+1.061*alfam*betam));
kg=nm*10^-7*(R*10^-3)*3.75*td/Mlm;
%------------------------------------------------------------------------
% Cálculo dos parâmetros hidrodinâmicos
%------------------------------------------------------------------------
Dex=26.7*10^-3;                      %Diâmetro externo do tubo
ntubos=49;
At=(pi*DI^2/4)-(pi*ntubos*Dex^2/4);   %Área da seção transversal do leito
Ar=dp^3*DENgl*(DENs-DENgl)*g/mi^2;    %Número de Arquimedes
esf=6/(dp*DENs*EAP);                  %Esfericidade
emf=0.586*esf^-0.72*Ar^-0.029*(DENgl/DENs)^0.021; % emf=0.586*(mi^2/(DENgl*g*(DENs-DENgl)*dp^3))^0.029*(DENgl/DENs)^0.021;
umf=(esf*dp)^2*g*(DENs-DENgl)*emf^3/(150*mi*(1-emf));             %Velocidade na mínima fluidização
uo=QF/At;                          %Velocidade superficial do gás
dbo=0.347*(Ao*(uo-umf))^0.4;          %Diâmetro mínimo das bolhas
dbm=0.652*(pi*DI^2*(uo-umf)/4)^0.4;     %Diâmetro máximo das bolhas
Zmf=W/(At*(1-emf)*(DENs-DENgl));      %Altura à mínima fluidização
db=dbm-(dbm-dbo)*exp(-0.3*Zmf/DI);      %Diâmetro da bolha
ub=uo-umf+0.711*(g*db)^0.5;          %Velocidade de ascenção das bolhas
delta=(uo-umf)/(ub-umf);              %Fração de bolhas no leito
ef=1-((1-delta)*(1-emf));             %Porosidade do leito fluidizado
Z=Zmf*(1-emf)/(1-ef);    %todo: fix   %Altura do leito fluidizado
Zproportion=Z/Zmf;                    %Relação da altura do leito sobre a mínima altura, deve ser 1,2
DeltaPf=150*mi*Z*uo*(1-ef)^2/(dp^2*ef^3)+1.75*Z*uo^2*DENgl*(1-ef)/(dp*ef^3);
AreaPerfurada=pi*DI^2/4*0.1;
Cd=0.85*(0.01/dorificio)^0.13;        %Coeficiente da equação bernoulli, 0.01 é a espessura do distribuidor
DeltaPDist=312.5*DENgl*(25*uo)^2/Cd^2;      %Perda de carga no distribuidor%
Ab=delta*At;                          %Área transversal da fase bolha
Qb=(uo-umf)*Ab;                       %Vazão volumétrica na fase bolha
Qe=QF-Qb;                             %Vazão volumétrica na fase emulsão
QeF=Qe;                               %Vazão volumétrica na fase emulsão na alimentação
V=At*Z;                               %Volume do leito fluidizado
Re=(DENgl*dp*umf)/mi;                 %Número de Reynolds
dstar=dp*(g*DENgl*(DENs-DENgl)/mi^2)^(1/3); %Diâmetro adimensional das partículas
ustar=1/(18/dstar^2+0.5909/dstar^0.5); %Velocidade terminal adimensional do gás
ut=ustar/(DENgl^2/(mi*(DENs-DENgl)*g))^(1/3); %Velocidade terminal do gás (m/s)
%Cd=18/Re^(3/5);
%ut=(4*(DENs-DENgl)*g*dp/(3*DENgl*Cd))^0.5;
Tres2=(At*Z)/QF;
UU=uo/umf;
Dleito=4*W/(pi*DI^2*Z);
Cc=0.04*DENs/134.4;                   %Concentração de cloreto de cobre no catalisador
%------------------------------------------------------------------------
% Cálculo dos coeficientes de transferência de massa
%------------------------------------------------------------------------
for i=1:6
    Kbc(i)=4.5*(umf/db)+5.85*((Di(i)^0.5*g^0.25)/db^1.25);              %Coef. de transf. de massa bolha/nuvem
end
for i=1:6
    Kce(i)=6.78*(emf*Di(i)*ub/db^3)^0.5;                                %Coef. de transf. de massa nuvem/emulsão
end
for i=1:6                          
    Kb(i)=1/Kbc(i)+(1/Kce(i));
    Kbe(i)=1/Kb(i);                                                     %Coef. de transf. de massa global
end
Hbc=4.5*(umf*DENgl*CPgl/db)+5.85*((kg*DENgl*CPgl)^0.5*g^0.25)/db^1.25;
Hce=6.78*((DENgl*CPgl*kg)^0.5)*(emf*ub/(db^3))^0.5;
Hbe=1/(1/Hbc+1/Hce);
for i=1:6
    teta(i)=Kbe(i)/(uo-umf);                                           %teta da equação de BM
end
teta1=teta(1);
teta2=teta(2);
teta3=teta(3);
teta4=teta(4);
teta5=teta(5);
teta6=teta(6);

%------------------------------------------------------------------------
% Cálculo da variável beta
%------------------------------------------------------------------------
beta=Hbe/(ub*DENgl*CPgl);

%------------------------------------------------------------------------
% Resolução do sistema de equações não lineares
%------------------------------------------------------------------------
x0=[N1e,N2e,N3e,N4e,N5e,N6e];
options=optimset('Display','iter','HessUpdate','bfgs');
[x,fval,exitflag,algorithm]=fsolve(@sistemaleito,x0,options,Qe,QF,N1F,N2F,N3F,N4F,N5F,N6F,Te,uo,umf,Ab,Z,V,teta1,teta2,teta3,teta4,teta5,teta6,ef,Kr,Ka,Cc,P,R,delta);
erro=(abs(x(1)-N1e)+abs(x(2)-N2e)+abs(x(3)-N3e)+abs(x(4)-N4e)+abs(x(5)-N5e)+abs(x(6)-N6e))/6;
if erro>0.0001
    N1e=x(1); N2e=x(2); N3e=x(3); N4e=x(4); N5e=x(5); N6e=x(6);
    conta_massa=conta_massa+1;
    erro1=erro;
    continue
end
if conta_massa>10
    N1e=x(1); N2e=x(2); N3e=x(3); N4e=x(4); N5e=x(5); N6e=x(6);
    break
end
erro2=erro;
N1e=x(1); N2e=x(2); N3e=x(3); N4e=x(4); N5e=x(5); N6e=x(6);
if N2e<0
    N2e=0;
end
%------------------------------------------------------------------------
% Recálculo da viscosidade, densidade e CP no leito e na saída para uso no BE
%------------------------------------------------------------------------
% Fator de compressibilidade---------------------------------------------
for i=1:6
    a(i)=0.42748*(R^2*Tc(i)^2.5/Pc(i));
    b(i)=0.08664*(R*Tc(i)/Pc(i));
end
bmF=y1F*b(1)+y2F*b(2)+y3F*b(3)+y4F*b(4)+y5F*b(5)+y6F*b(6);
bme=y1e*b(1)+y2e*b(2)+y3e*b(3)+y4e*b(4)+y5e*b(5)+y6e*b(6);
amF=(y1F*(a(1)^0.5)+y2F*(a(2)^0.5)+y3F*(a(3)^0.5)+y4F*(a(4)^0.5)+y5F*(a(5)^0.5)+y6F*(a(6)^0.5))^2;
ame=(y1e*(a(1)^0.5)+y2e*(a(2)^0.5)+y3e*(a(3)^0.5)+y4e*(a(4)^0.5)+y5e*(a(5)^0.5)+y6e*(a(6)^0.5))^2;
VidF=R*TMF/PMF;
Vide=R*Te/P;
VcalcF=@(VF)(PMF+amF/(TMF^0.5*(VF+bmF)*VF))*(VF-bmF)-R*TMF;
Vcalce=@(Ve)(P+ame/(Te^0.5*(Ve^2+bme*Ve)))*(Ve-bme)-R*Te;
VrF=fzero(VcalcF,VidF);
Vre=fzero(Vcalce,Vide);
ZF=VrF/VidF;
Ze=Vre/Vide;
% Viscosidade no leito---------------------------------------------------
Mm=y1e*M(1)+y2e*M(2)+y3e*M(3)+y4e*M(4)+y5e*M(5)+y6e*M(6);
Tcm=y1e*Tc(1)+y2e*Tc(2)+y3e*Tc(3)+y4e*Tc(4)+y5e*Tc(5)+y6e*Tc(6);
Zcm=y1e*Zc(1)+y2e*Zc(2)+y3e*Zc(3)+y4e*Zc(4)+y5e*Zc(5)+y6e*Zc(6);
Vcm=(y1e*Vc(1)+y2e*Vc(2)+y3e*Vc(3)+y4e*Vc(4)+y5e*Vc(5)+y6e*Vc(6))*1000;
Pcm=R*10^-2*Tcm*(Zcm/Vcm);
Trm=Te/Tcm;
Prm=P*10^-5/Pcm;
zeta=0.176*(Tcm/(Mm^3*Pcm^4))^(1/6);
for i=1:6
    Tr(i)=Te/Tc(i);
end
for i=1:6
    Dpr(i)=(52.46*Dp(i)^2*Pc(i))/Tc(i)^2;
end
for i=1:6
    if Dpr(i)>=0 & Dpr(i)<0.022
        Fpo(i)=1;
    elseif Dpr(i)>=0.022 & Dpr(i)<0.075
        Fpo(i)=1+30.55*(0.292-Zc(i))^1.72;
    elseif Dpr(i)>=0.075
        Fpo(i)=1+30.55*(0.292-Zc(i))^1.72*abs((0.96+0.1*(Tr(i)-0.7)));
    end
end
Fpom=y1e*Fpo(1)+y2e*Fpo(2)+y3e*Fpo(3)+y4e*Fpo(4)+y5e*Fpo(5)+y6e*Fpo(6);
Z1=((0.807*Trm^0.618-0.357*exp(-0.449*Trm)+0.340*exp(-4.058*Trm)+0.018))*Fpom;
Visco=Z1/zeta;
mi=Visco*10^-7;
% Capacidades caloríficas na entrada, leito e saída----------------------
for i=1:6
    CPF(i)=(CP(i,1)+CP(i,2)*10^-3*TMF+CP(i,3)*10^-5*TMF^2+CP(i,4)*10^-8*TMF^3+CP(i,5)*10^-11*TMF^4)*R;
end
for i=1:6
    CPo(i)=(CP(i,1)+CP(i,2)*10^-3*512.15+CP(i,3)*10^-5*512.15^2+CP(i,4)*10^-8*512.15^3+CP(i,5)*10^-11*512.15^4)*R;
end
for i=1:6
    CPl(i)=(CP(i,1)+CP(i,2)*10^-3*Te+CP(i,3)*10^-5*Te^2+CP(i,4)*10^-8*Te^3+CP(i,5)*10^-11*Te^4)*R;
end
CPFar=(3.355+0.575*10^-3*TMF-0.016*10^5*TMF^-2)*R;
CPgF=(CPF(1)*y1F+CPF(2)*y2F+CPFar*yar);
CPgo=(CPo(1)*y1o+CPo(2)*y2o+CPo(3)*y3o+CPo(4)*y4o+CPo(5)*y5o+CPo(6)*y6o);
CPgl=(CPl(1)*y1e+CPl(2)*y2e+CPl(3)*y3e+CPl(4)*y4e+CPl(5)*y5e+CPl(6)*y6e);
% BM na fase bolha-------------------------------------------------------
N1b=((N1F/QF-N1e/Qe)*exp(-teta1*Z)+N1e/Qe)*Qb;
N2b=((N2F/QF-N2e/Qe)*exp(-teta2*Z)+N2e/Qe)*Qb;
N3b=((N3F/QF-N3e/Qe)*exp(-teta3*Z)+N3e/Qe)*Qb;
N4b=((N4F/QF-N4e/Qe)*exp(-teta4*Z)+N4e/Qe)*Qb;
N5b=((N5F/QF-N5e/Qe)*exp(-teta5*Z)+N5e/Qe)*Qb;
N6b=((N6F/QF-N6e/Qe)*exp(-teta6*Z)+N6e/Qe)*Qb;
% Fluxo na saída---------------------------------------------------------
Nout(1)=N1b+N1e;     Nout(2)=N2b+N2e;     Nout(3)=N3b+N3e;
Nout(4)=N4b+N4e;     Nout(5)=N5b+N5e;     Nout(6)=N6b+N6e;
% Densidade--------------------------------------------------------------
MMe=M(1)*y1e+M(2)*y2e+M(3)*y3e+M(4)*y4e+M(5)*y5e+M(6)*y6e;
MMo=M(1)*y1o+M(2)*y2o+M(3)*y3o+M(4)*y4o+M(5)*y5o+M(6)*y6o;
DENgl=P*MMe/(Ze*R*Te);
DENgo=P*MMo/(Ze*R*Te);
%------------------------------------------------------------------------
% Cálculo do calor de reação
%------------------------------------------------------------------------
Da=(CP(3,1)+CP(5,1)-2*CP(2,1)-0.5*CP(4,1)-CP(1,1))*R;
Db=(CP(3,2)+CP(5,2)-2*CP(2,2)-0.5*CP(4,2)-CP(1,2))*R*10^-3;
Dc=(CP(3,3)+CP(5,3)-2*CP(2,3)-0.5*CP(4,3)-CP(1,3))*R*10^-5;
Dd=(CP(3,4)+CP(5,4)-2*CP(2,4)-0.5*CP(4,4)-CP(1,4))*R*10^-8;
De=(CP(3,5)+CP(5,5)-2*CP(2,5)-0.5*CP(4,5)-CP(1,5))*R*10^-11;
DH=(-2.39*10^8+Da*(Te-Tref)+Db*0.5*(Te^2-Tref^2)+Dc*0.33*(Te^3-Tref^3)+Dd*0.25*(Te^4-Tref^4)+De*0.2*(Te^5-Tref^5));
r1=(Kr*Ka*Cc*N1e*P/(N1e+N2e+N3e+N4e+N5e+N6e))/(R*Te+(Ka*N1e*P/(N1e+N2e+N3e+N4e+N5e+N6e)));
%------------------------------------------------------------------------
% Cálculo do coeficiente global de troca térmica
%------------------------------------------------------------------------
Dex=26.7*10^-3;    Din=20.96*10^-3;
L=(DI*1.5+0.865)*ntubos; % 0.865 vem da curva no fim do tubo
Kt=70;
Cr=0.05;            Pr=0.902;           magua=43;
Ain=(pi*Din*L);
DENagua=992.2;
uagua=(magua/DENagua)/(ntubos*Din^2*pi/4);
Aw=(pi*Dex*L);
Rew=(uagua*DENagua*Din)/(547*10^-6);
Nu=0.023*(Rew^0.8)*(Pr^0.33);
hi=(Nu*0.647)/Din;
hw=0.88*kg/dp*(DENgl*(DENs-DENgl)*g*dp^3/mi^2)^0.213; %hw=0.01844*Cr*(kg/dp)*(1-ef)*(DENgl*CPgl/kg)^0.43*((dp*DENgl*uo)/mi)^0.23*(417/CPgl)^0.8*(DENs/DENgl)^0.66;
rt=(1/(hi*Ain))+(log(Dex/Din)/(2*pi*Kt*L))+(1/(hw*Aw));
U=1/(rt*Aw);
%------------------------------------------------------------------------
% BE na fase emulsão
%------------------------------------------------------------------------
Qo=QF;
Qeo=Qo-Qb;
CalorAgua=U*Aw*(Tw-Te);
DeltaTagua=-CalorAgua/(magua*4180);
Energia=QeF*DENgF*(CPgF/MMF)*(TMF-Tref)-Qeo*DENgo*(CPgo/MMo)*(Te-Tref)+ub*DENgl*Ab*(CPgl/MMe)*(TMF-Te)*(1-(exp(-beta*Z)))+At*Z*(1-delta)*(1-ef)*(-1)*DH*r1+U*Aw*(Tw-Te);
if abs(Energia)>1
    Te=(QeF*DENgF*(CPgF/MMF)*(TMF-Tref)-Qeo*DENgo*(CPgo/MMo)*(-Tref)+ub*DENgl*Ab*(CPgl/MMe)*TMF*(1-(exp(-beta*Z)))+At*Z*(1-delta)*(1-ef)*(-1)*DH*r1+U*Aw*Tw)/(Qeo*DENgo*(CPgo/MMo)+ub*DENgl*Ab*(CPgl/MMe)*(1-(exp(-beta*Z)))+U*Aw);
    loop_energia=1;
end
if abs(Energia)<1
    loop_Energia=0;
end
if conta_energia==10
    loop_energia=0;
end
conta_energia=conta_energia+1;
end                            %Fechamento do while do início do programa
%------------------------------------------------------------------------
% Cálculo das conversões na região de leito
%------------------------------------------------------------------------
Conv.etilenoleito=((N1F-Nout(1))/N1F)*100;
Conv.HClleito=((N2F-Nout(2))/N2F)*100;
Conv.O2leito=((N4F-Nout(4))/N4F)*100;
fprintf('Massa na saída do leito:             %6.5e\n', Nout(1)*M(1)+Nout(2)*M(2)+Nout(3)*M(3)+Nout(4)*M(4)+Nout(5)*M(5)+Nout(6)*M(6));
fprintf('Massa de N2 na entrada do leito:       %6.3e\n', mN2entrada);
fprintf('Massa de N2 na saída do leito:       %6.3e\n', Nout(6)*M(6));
%------------------------------------------------------------------------
% Cálculo da região de freeboard
%------------------------------------------------------------------------
N6=Nout(6);

% cálculo do uo no freeboard
Q1FreeBoard=Nout(1)*R*Te/P;
Q2FreeBoard=Nout(2)*R*Te/P;
Q3FreeBoard=Nout(3)*R*Te/P;
Q4FreeBoard=Nout(4)*R*Te/P;
Q5FreeBoard=Nout(5)*R*Te/P;
Q6FreeBoard=Nout(6)*R*Te/P;
QFFreeBoard=(Q1FreeBoard+Q2FreeBoard+Q3FreeBoard+Q4FreeBoard+Q5FreeBoard+Q6FreeBoard);
uoFreeBoard=QFFreeBoard/(pi*DI^2/4);
fprintf('uoFreeBoard');
disp(uoFreeBoard);

a=0.75/uoFreeBoard; % 0.75 é a 'Constante de decaimento' (Kunii, D.; Levenspiel, O., 1990)
A=(pi*DI^2)/4;
s1=Nout(1);  s2=Nout(2);  s3=Nout(3);  s4=Nout(4);  s5=Nout(5);  s6=Nout(6);
% Resolução do sistema de equações diferenciais--------------------------
options=odeset('RelTol',1e-5,'AbsTol',[1e-5]);
Zfinal=7.5-Z;
altfreeboard=0:0.05:Zfinal;
disp(Z);
disp(Zfinal);
[Zf,N]=ode15s(@sistemafreeboard,[altfreeboard],[s1,s2,s3,s4,s5],options,ef,A,a,Kr,Ka,Cc,P,R,Te,s6);
% Z,N,ef,A,a,Kr,Ka,Cc,P,R,T,N6
%[x,fval,exitflag,algorithm]=fsolve(@sistemaleito,x0,options,Qe,QF,N1F,N2F,N3F,N4F,N5F,N6F,Te,ub,uo,umf,Ab,Z,V,teta1,teta2,teta3,teta4,teta5,teta6,ef,Kr,Ka,Cc,P,R,DENgo,CPgo,beta,delta);
%x,Qe,QF,N1F,N2F,N3F,N4F,N5F,N6F,Te,ub,uo,umf,Ab,Z,V,teta1,teta2,teta3,teta4,teta5,teta6,ef,Kr,Ka,Cc,P,R,DENgo,Cpgo,beta,delta
linha=length(N);
%disp(N);
matN=1;
for i=1:linha
    matN=i;
    if N(i,2)<0.00001
        break
    end
end
if matN < 1
    matN = 1;
end
N1saida=N(matN,1);
N2saida=N(matN,2);
N3saida=N(matN,3);
N4saida=N(matN,4);
N5saida=N(matN,5);
N6saida=s6; 
%  Conversões na região de freeboard-------------------------------------
Conv.etilenofreeboard=((N1F-N1saida)/N1F)*100;
Conv.HClfreeboard=((N2F-N2saida)/N2F)*100;
Conv.O2freeboard=((N4F-N4saida)/N4F)*100;
QSaida=(N1saida+N2saida+N3saida+N4saida+N5saida+N6saida)*R*Te/P;
DENsaida=(P/(R*Te))*(N1saida*M(1)+N2saida*M(2)+N3saida*M(3)+N4saida*M(4)+N5saida*M(5)+N6saida*M(6))/(N1saida+N2saida+N3saida+N4saida+N5saida+N6saida);
%
%---------------------------- FIM DO PROGRAMA ---------------------------
%
%------------------------------------------------------------------------
% Impressão dos resultados
%------------------------------------------------------------------------
fprintf('%s\n','-------------------------------------------------------')
fprintf('%s\n','Resultados para o leito')
fprintf('%s\n','-------------------------------------------------------')
fprintf('Taxa de conversao do etileno:  %6.2f\n',Conv.etilenoleito)
fprintf('Taxa de conversao do HCL:      %6.2f\n',Conv.HClleito)
fprintf('Taxa de conversao do Oxigênio: %6.2f\n',Conv.O2leito)
fprintf('%s\n','-------------------------------------------------------')
fprintf('%s\n','Resultados para o leito + freeboard')
fprintf('%s\n','-------------------------------------------------------')
fprintf('Taxa de conversao do etileno:    %6.2f\n',Conv.etilenofreeboard)
fprintf('Taxa de conversao do HCl:        %6.2f\n',Conv.HClfreeboard)
fprintf('Taxa de conversao do Oxigênio:   %6.2f\n',Conv.O2freeboard)
fprintf('%s\n','-------------------------------------------------------')
fprintf('%s\n','Resultados parametros hidrodinamicos e termodinamicos do leito')
fprintf('%s\n','-------------------------------------------------------')
fprintf('Porosidade do leito expandido:                 %6.3f\n',ef)
fprintf('Porosidade do leito na minima fluidizaçao:     %6.3f\n',emf)
fprintf('Altura do leito expandido:                     %6.3f\n',Z)
fprintf('Altura do leito na minima fluidizacao          %6.3f\n',Zmf)
fprintf('Relação da altura do leito                     %6.3f\n',Z/Zmf)
fprintf('Velocidade de minima fluidiacao:               %6.3e\n',umf)
fprintf('Velocidade no Leito (velocidade superficial):  %6.3f\n',uo)
fprintf('Velocidade terminal das partículas:            %6.3e\n',ut)
fprintf('Velocidade de ascençao da bolha:               %6.3f\n',ub)
fprintf('Temperatura no leito:                          %6.2f\n',Te)
fprintf('Perda de Pressão no Leito:                     %6.3f\n',DeltaPf)
fprintf('Perda de Pressão no Distribuidor:              %6.3f\n',DeltaPDist)
fprintf('N° de orificios:                               %6.3f\n',norificios)
fprintf('Diametro do orificio:                          %6.3f\n',dorificio)
fprintf('Pitch:                                         %6.3f\n',pitch)
fprintf('Diametro mínimo da bolha:                      %6.3f\n',dbo)
fprintf('Diametro medio da bolha:                       %6.3f\n',db)
fprintf('Diametro máximo da bolha:                      %6.3f\n',dbm)
fprintf('Viscosidade dos gases no leito:                %6.3e\n',mi)
fprintf('Condutividade termica dos gases no leito:      %6.3e\n',kg)
fprintf('Calor especifico da mistura entrada reator:    %6.3e\n',CPgF)
fprintf('Calor especifico da mistura saida reator:      %6.3e\n',CPgo)
fprintf('Calor especifico da mistura no leito:          %6.3e\n',CPgl)
fprintf('Velocidade da água de resfriamento:            %6.3e\n',uagua)
fprintf('DeltaT da água de resfriamento:                %6.3e\n',DeltaTagua)
fprintf('Vazão volumétrica de produto:                  %6.3e\n',QSaida)
fprintf('Vazão volumétrica da entrada:                  %6.3e\n',Qentrada)
fprintf('Viscosidade dos gases:                         %6.8f\n',mi)
fprintf('Densidade dos gases na entrada do reator:      %6.3f\n',DENentrada)
fprintf('Densidade dos gases na saida do reator:        %6.3f\n',DENsaida)
fprintf('Densidade dos gases no leito do reator:        %6.3f\n',DENgl)
fprintf('Área externa dos tubos do trocado:             %6.3f\n',Aw)
fprintf('Coeficiente de troca térmica:                  %6.3f\n',U)