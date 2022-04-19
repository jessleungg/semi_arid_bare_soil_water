function unconfinedaquifer

%this function is to stimulate the movement of water and vapour through a vertical 
%column in a semi-arid environment with an initial wetting front and a 
%horizontal distance of L. A zero flux boundary has been used and times(t) 
%have been analysed with finite difference. 


%Van genuchten paramaters for the soil in SI units where applicable, if no
%units stated variable is unitless
thS=0.469;
thR=0.19;
Ks=303/100/24/60^2; %(ms-1)
alpha=0.005*100; %(m-1)
n=7.09;
m=1-1/n;
eta=0.5;
rhow = 997; %density of water (kgm-3)


%times of interest
t =[0 10 20 30 40 50]*60*100000; %(s)
L = 10; %(m)


%Initial and boundary conditions determined from a stated RH  
RH=[0.95]'; %relative humidity
T = 290.75; %temperature (K)
Rw = 461.4; %individual gas constant for water (Jkg-1K-1)
g = 9.81; %gravity (ms-2)
psi0=log(RH)*Rw*T/g; %pressure head (m)


%spatial grid
N = 100;
xB=[0 logspace(-8,0,N)]'*L; %x boundaries - creates a log spaced column of N+1 values from 0-L
x = (xB(1:N,1) + xB(2:N+1,1))/2; %x midpoints (m)


SeI=0.99; %Initial effective saturation
psiI = -(SeI.^(-1/m)-1).^(1/n)/alpha; %initial pressure head (m)
psiI_vec = zeros(size(x)) + psiI; %initial pressure head vector


%Initial conditions for cummulative mass flux
F0I=0;
FNI=0;
psiI_vec=[F0I;psiI_vec;FNI]; %Add these to psiI_vec


JPat=spdiags(ones(N+2,3),[-1 0 1],N+2,N+2); %initialising results matrix
options=odeset('JPattern',JPat);
[t, psiFD] = ode15s(@MYodefun,t,psiI_vec,options,psi0,x,xB,alpha,n,m,Ks,eta,thS,thR,rhow,RH); 


for k = 1:numel(t) %extracting mass flux of vapour and specific moisture capacity from ODE
[~,Fv(:,k),dmdpsi(:,k),F0(k),FN(k),Mw(k)] = MYodefun(t(k),psiFD(k,:)',psi0,x,xB,alpha,n,m,Ks,eta,thS,thR,rhow,RH) ;
end


%converting pressure head to moisture content for FD
psiFD(:,[1 end])=[]; %Take out cummulative mass fluxes
%Apply Van Genuchten soil moisture characteristics
Se=(1+abs(alpha*psiFD).^n).^-m; %effective saturation
Ss=0; %Specific Storage Coefficient (m-1)
dSedpsi=m./(1-m).*Se.*(Se.^(1/m)-1)./psiFD; %(m-1)
ind=psiFD>=0;
Se(ind)=1;
dSedpsi(ind)=0;
%Apply Webb extension for oven drying
[Se,dSedpsi]=WebbExtension(psiFD,Se,dSedpsi,thR,thS,alpha,n,m);
thetaFD=Se.*(thS-thR)+thR; %moisture content - finite difference


%plotting figures
figure(1)
clf
subplot(2,1,1)
hold off
plot(x,Se,'.-')
hold on
set(gca,'ColorOrderIndex',1)
xlabel('Distance (m)')
ylabel('Effective saturation (-)')
set(gca,'xscale','log')
subplot(2,1,2)
plot(t,Mw-Mw(1),'k',t,F0-FN,'ro')
xlabel('Time (s)')
ylabel('Cummulative water mass')

%**************************************************************************

function [dpsidt,Fv,dmdpsi,F0,FN,Mw] = MYodefun(t,psi,psi0,x,xB,alpha,n,m,Ks,eta,thS,thR,rhow,RH) 

%Extract cummulative fluxes from solution vector
F0=psi(1);
FN=psi(end);
psi(1)=[];
psi(end)=[];
x = [xB(1);x];
psi = [psi0;psi];


%Apply Van Genuchten soil moisture characteristics
Se=(1+abs(alpha*psi).^n).^-m; %effective saturation
Ss=0; %Specific Storage Coefficient (m-1)
dSedpsi=m./(1-m).*Se.*(Se.^(1/m)-1)./psi; %(m-1)
ind=psi>=0;
Se(ind)=1;
dSedpsi(ind)=0; %(m-1)

K=Ks*Se.^eta.*(1-(1-Se.^(1/m)).^m).^2; %hydraulic conductivity (ms-1)
K(Se<0)=0;


%Apply Webb extension for oven drying
[Se,dSedpsi]=WebbExtension(psi,Se,dSedpsi,thR,thS,alpha,n,m);
theta=Se.*(thS-thR)+thR; %moisture content

KB =interp1q(x,K,xB(1:end-1)); %hydraulic conductivity
KB(1) = K(1); %boundary condition


%Darcy's law
dpsidx = diff(psi,1,1)./diff(x,1,1);
q = -KB .* dpsidx; %Darcy flux (ms-1)
q = [q;0]; %zero flux boundary

%removing boundaries
psi(1) = []; 
theta(1) = [];
x(1)=[];


%vapour component
[Fv,Cv,dmdpsi] = VapourStuffFun(psi,theta,thS,RH,x,xB,thR,dSedpsi,rhow); %include vapour component
q = q + Fv/rhow; %adding to vapour to Darcy's flux (ms-1)
C = dmdpsi/rhow; %specific moisture capacity (m-1)
mw=rhow*theta+(thS-theta).*Cv; %mass of H2O per unit volume of rock (kgm-3)
Mw=sum(mw.*diff(xB)); %cummulative water mass


%calculate dpsidt with vapour
dqdx = diff(q,1,1)./diff(xB,1,1); %(s-1)
dthetadt = -dqdx; %mass conservation law (s-1)
dpsidt = dthetadt./C; %applying product rule (ms-1)


%Determine derivatives of cummulative mass fluxes
dF0dt=rhow*q(1);
dFNdt=rhow*q(end);
%Add these to the time-derivative vector
dpsidt=[dF0dt;dpsidt;dFNdt];

%**************************************************************************

function [Se,dSedpsi]=WebbExtension(psi,Se,dSedpsi,thR,thS,alpha,n,m)

%Determine dSwdpsi
Swr=thR/thS; %residual water saturation
Sw=(1-Swr)*Se+Swr; %water saturation
dSwdpsi=(1-Swr)*dSedpsi; %(m-1)
%Apply Web drying
psiD=-1e5; %(m)
psiDD=-alpha*psiD;
%Swm - the water saturation below which the webb modle applies
[Swm,psiMD]=WebbSem(m,n,Swr,psiDD); %to when Webb is applied
%calculate Sw if the webb extension is applied 
SwWEB=Swm*log(psi/psiD)/log(psiMD/psiDD); %Webb water saturation
dSwdpsiWEB=Swm./psi./log(psiMD/psiDD); %Webb dSwdpsi
%acts like an if function
%if sw<swm then Sw for those indexed values will be made equal to swWebb
%else Sw above those values retains the original van genuchten value0
ind=Sw<=Swm;
Sw(ind)=SwWEB(ind);
dSwdpsi(ind)=dSwdpsiWEB(ind);
%redefine dSedpsi and Se - converting water saturation to effective
%saturation (-ve)
dSedpsi=dSwdpsi/(1-Swr);
Se=(Sw-Swr)/(1-Swr);

%**************************************************************************

function [Sm,Pm]=WebbSem(m,n,Sc,Pd)
z=Sc./(1-Sc).*Pd.^(m.*n); %lambert variable
%lambert W function
L1=log(z); %lambert variable
L2=log(L1); %lambert variable
W=real(L1-L2+L2./L1);
Sem=Sc./(1-Sc)./W; %first intial guess for sem - critical efective saturation 
%iterative process
for k=1:10
   Pm=(Sem.^(-1./m)-1).^(1./n); %capillary pressure at critical effective saturation (kgm-1s-2)
   Sm=(1-Sc).*Sem+Sc; %aqueous saturation matching point 
   Sem=Sm./(1-Sc)./n./m./(Sem.^(1./m)-1)./log(Pm./Pd); %critical effetive saturation
end

%**************************************************************************

function [Fv,Cv,dmdpsi]=VapourStuffFun(psi,theta,thS,RH,x,xB,thR,dSedpsi,rhow)

%initialising parameters
T = 290.75; %temperature (K)
Rw = 461.4; %individual gas constant for water (Jkg-1K-1)
P0 = 0; %atmospheric pressure (Nm-2)
g = 9.81; %gravity (ms-2)
DEa = 2.5e-5; %effective diffusion of vapour in the air (m2s-1)
tort=1; %tortuosity factor

%Ficks first Law - Mass flux of Vapour
Sv = 1 - theta/thS; %volume fraction of pore space occupied by vapour
Sw = 1 - Sv; %volume fraction of pore space occupied by water
phi = thS;

Pvs = 610.78 * exp(17.27 * (T - 273.15)/(T - 35.85)); % water pressure at saturation - tetens equation (kgm-1K-2)
Pv = Pvs * exp(psi*g/Rw/T); %vapour pressure (kgm-1K-2)

De = tort*phi.*Sv*DEa; %effective diffusion coefficient (m2s-1)
De = interp1q(x,De,xB(1:end-1)); %boundary coefficients
De(1)=De(2);

Cv = Pv/Rw/T; %mass concentration of water vapor (kgm-3)
Cv0=RH*Pvs/Rw/T; %boundary condition

dCvdx = diff([Cv0;Cv],1,1)./diff([xB(1);x],1,1);
Fv =[-De .*dCvdx;0]; %Ficks vapour law (kgm-2s-1)


%specific capacity
Pw0 = 1e5; %boundary water pressure (kgm-1s-2)
Pw = Pw0 + rhow * g * psi; %water pressure
Cw = rhow; %mass concentration of water (kgm-3)
cV = 1/rhow/Rw/T; %compressibility of vapour (ms2kg-1)

cR = 0; %compressibility of rock (ms2kg-1)
cW = 4.4e-10; %compressibility of water (ms2kg-1)
dSedpsi(1) = []; %remove boundary condition

dmdpsi = ((Cw - Cv)*(thS - thR) .* dSedpsi) + Cw*g*thS.*(Cw.*Sw.*(cW + cR) + Cv.*Sv.*(cV + cR)); %specific capacity (kgm-4)


