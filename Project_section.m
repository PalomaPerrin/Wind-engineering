%%Wind Eng Lab - Bridge Aeroelasticity
clc
clear 
close all


%Constants
rho = 1.25; % air density [kg/m3]   
B   = 35.6; % deck chord [m]
L   = 1;    % [m] --> unit length

% Structural characteristics 
my = 23160;       % [kg/m]
mz = my;          % 
J  = 2.77e6;      % [kg m^2 / m] 

zeta  = 4e-3;     % damping ratio:  zeta = r / (2 m omega)

fy = 0.05;        %natural frequencies [Hz]
fz = 0.0884;
ft = 0.259;

ky = my*(2*pi*fy)^2;  ry = zeta*2*my*(2*pi*fy);
kz = mz*(2*pi*fz)^2;  rz = zeta*2*mz*(2*pi*fz);
kt = J*(2*pi*ft)^2;   rt = zeta*2*J*(2*pi*ft);

%Structural Matrices
M_stru = diag([my,mz,J]);
R_stru = diag([ry,rz,rt]);
K_stru = diag([ky,kz,kt]);

% Aero Coefficents
load('IBB_polari'); %alpha is in [deg] !
Cd = DLM(:,1);    %Tutte le righe della prima colonna del file DLM: ...
                  %Contiene i valori di Cd, Cl, Cm per diversi alpha 
Cl = DLM(:,2);
Cm = DLM(:,3);
Kd = KDLM(:,1);   %File KDLM: contiene i valori delle pendenze ... slope of the coefs
                  %Di Cd, Cl, Cm per diversi alpha 
Kl = KDLM(:,2);
Km = KDLM(:,3);

figure (1)
subplot(1,2,1)
plot(alpha,Cd,'-or','LineWidth',2)
hold on
plot(alpha,Cl,'-og','LineWidth',2)
plot(alpha,Cm,'-ob','LineWidth',2)
legend('CD','CL','CM');
grid on; title('D,L,M')   
xlabel('Angle of attack \alpha [deg]','FontSize',10)
ylabel('[Coefficents]','FontSize',10)

subplot(1,2,2)
plot(alpha,Kd,'-or','LineWidth',2)
hold on
plot(alpha,Kl,'-og','LineWidth',2)
plot(alpha,Km,'-ob','LineWidth',2)
legend('KD','KL','KM');
grid on; title('D,L,M derivatives')   
xlabel('Angle of attack \alpha [deg]','FontSize',10)
ylabel('[Slopes]','FontSize',10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATIC SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
U_min = 0;  U_max = 80;  nvel = 80;
U = linspace(U_min,U_max,nvel);  %wind speed vector 
theta0 = [-10:0.1:10];    %Prendiamo un range di rotazioni statiche theta0

%Simplified static solution
CM0 = interp1(alpha,Cm,0);   %Non avendo alpha = 0 nell'array di valori faccio...
                             %L'interpolazione coi valori vicini per trovare Cm0
KM0 = interp1(alpha,Km,0);

CM_lin = CM0 + KM0*theta0*pi/180;  %Cm linearizzato

figure (2)
plot(theta0, CM_lin,'-r','LineWidth',2,'DisplayName','C_{M0} + K_{M0} \theta_{0}')
hold on
grid on
plot(alpha,Cm,'-ob','LineWidth',1)
legend('C_{M0} + K_{M0}*\theta_0', 'Cm(\theta)','Location','northwest')
xlabel('Angle of attack \theta_{0} [deg]','FontSize',10)
ylabel('Torsional coefficient','FontSize',10)

%
%Linearized solution:          0.5*rho*(U.^2)*(B^2)*CM0 
%                  theta_0 = ----------------------------
%                            kt -0.5*rho*(U.^2)*(B^2)*KM0

th_simple = (0.5*rho*(U.^2)*(B^2)*CM0)./(kt -0.5*rho*(U.^2)*(B^2)*KM0);

%Nonlinear solution
sol_old = 0;
t_statico = zeros(1,length(U));
errore_st = zeros(1,length(U));
hh = waitbar(0,'Computing static displacemnts (section model)...');

steps = length(U);
for ii = 1:length(U)
    u = U(ii);
        
    % NB : THETA and ALPHA are IN DEG, torsional stiffness is in Nm/rad
    [equilibrio] = fsolve(@(theta) statica(K_stru,theta,rho,u,B,L,alpha,Cm),sol_old);
    t_statico(ii) = equilibrio; %Theta static, in deg
    sol_old = equilibrio;       
    waitbar(ii/steps)
end
close(hh)

figure ('Name','Static displacements and angles ')

subplot(3,1,3)
plot(U,t_statico,'b','LineWidth',2);
hold on;
grid on;
subplot(3,1,3)
plot(U,th_simple*180/pi,'r');  %Lo passiamo da radianti a gradi
ylabel('\theta [deg]')
xlabel('Velocity U [m/s]')
legend('Non Linear','Linear','Location','northwest')
title('Static torsional displacement')

%Y and Z displacements: le rispettive matrici K sono per unità di lunghezzza!...
%Vanno ricavati moltiplicandole per B = 35.6
y_statico = 0.5*rho*(U.^2)*(B^2).*interp1(alpha,Cd,t_statico)/(B*ky);
z_statico = 0.5*rho*(U.^2)*(B^2).*interp1(alpha,Cl,t_statico)/(B*kz);

subplot(3,1,1)
plot(U,y_statico, 'r', 'LineWidth',2)
xlabel('Velocity U [m/s)')
ylabel('y [m]')
title('Static Y-axis displacement')
grid on;

subplot(3,1,2)
plot(U,z_statico, 'g', 'LineWidth',2)
xlabel('Velocity U [m/s)')
ylabel('z [m]')
title('Static Z-axis displacement')
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FLUTTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1y = 0*B; 
B1z = 0*B;
B1t = 0.5*B;

smo_y = zeros(1,length(U)); %damping ç 
smo_z = zeros(1,length(U));
smo_t = zeros(1,length(U));

freq_y = zeros(1,length(U)); %frequency [Hz]
freq_z = zeros(1,length(U));
freq_t = zeros(1,length(U));

autov = zeros(6,length(U));

hh = waitbar(0,'Calcolo velocità critica di flutter (SEZIONALE)...');
steps = length(U);

for ii = 1:length(U)
     u = U(ii); 
     qBL = 0.5*rho*(u^2)*B*L ;   
     t_0 = t_statico(ii);     %Estrapolo le condizioni statiche per ogni velocità
                              %theta0 -->  i Cl/k/m0;  le pend.ze Kd/l/m0
     Cd_eq = interp1(alpha,Cd,t_0);
     Cl_eq = interp1(alpha,Cl,t_0);
     Cm_eq = interp1(alpha,Cm,t_0);
     Kd_eq = interp1(alpha,Kd,t_0);
     Kl_eq = interp1(alpha,Kl,t_0);
     Km_eq = interp1(alpha,Km,t_0);

     %Aero Matrices
     %M_aero = zeros(3,3);   %La creiamo matrice di zero perchè non influenzata
     %R_aero = qBL.*[2*Cd_eq, (Kd_eq - Cl_eq), 0; ...
     %              2*Cl_eq, (Kl_eq + Cd_eq),  0;
      %             2*Cm_eq*B,    Km_eq*B,         Km_eq*B*B1t   ];
               
     %K_aero = -qBL.*[0, 0, Kd_eq; ...
       %             0, 0, Kl_eq;
        %            0, 0, Km_eq*B];
    M_aero = zeros(3,3);

    R_aero = [0.5*rho*u*B*2*Cd_eq  0.5*rho*u*B*(Kd_eq-Cl_eq)  0;
                   0.5*rho*u*B*2*Cl_eq  0.5*rho*u*B*(Kl_eq+Cd_eq)  0;
                   0.5*rho*u*B*B*2*Cm_eq 0.5*rho*u*B*B*Km_eq 0.5*rho*u*B*B*Km_eq*B1t];
             
    K_aero = [0 0 -0.5*rho*u^2*B*Kd_eq;
                   0 0 -0.5*rho*u^2*B*Kl_eq;
                   0 0 -0.5*rho*u^2*B^2*Km_eq];
     %Total matrices
     M = M_stru + M_aero;  %Le matrici 'stru' sono per unità di lungh.za!     
     R = R_stru + R_aero;  %Vanno moltiplicate per B
     K = K_stru + K_aero;
     
     % polyeig: eigenvalue/eigenvector solution
     [autovet,autoval] = polyeig(K,R,M); %Svolge ([M]lambda^2+([R]+R_aero)lambda...
                                         %         +([K]+K_aero)phi = 0
     
     %Riordino gli autovalori confrontando ognuno con quello al passo
     %precedente (in termini di modulo)
     if ii == 1
        f_old = [fy fz ft];
        smo_old = [zeta zeta zeta];
        lambda_old(1:3) = -smo_old.* (2*pi*f_old) + 1j* (2*pi*f_old) .* sqrt(1- smo_old.^2);
        lambda_old(4:6) = -smo_old.* (2*pi*f_old) - 1j* (2*pi*f_old) .* sqrt(1- smo_old.^2);
        %lambda_old è l'autovalore ottenuto dal sistema smorzato: 
        %lambda_old = - ç*omega_n +- j*omega_n*radice(1 - ç^2) 
        %           = - alpha +- j*omega_d
     end
     indici = 1:6;
     for jj = 1:6
         [val, dove] = min(abs(autoval(indici) - lambda_old(jj)));
         %'dove' è l'indice del valore più piccolo tra i minimi calcolati
         ind(jj) = indici(dove);
         indici = setdiff(indici,ind(jj));
     end
     autoval = autoval(ind);
     autovet = autovet(:,ind); 
     
     %damping:      Re(lambda)    -ç*(2*pi*f)
     %        ç = - --------- = - ----------
     %              |lambda|        2*pi*f
     smo_y(ii)= -real(autoval(1))/norm(autoval(1));
     smo_z(ii)= -real(autoval(2))/norm(autoval(2));
     smo_t(ii)= -real(autoval(3))/norm(autoval(3));
     %damped frequency:    |Im(lambda)|  2*pi*f*sqrt(1-ç^2)
     %                f_d = ---------- = ------------------
     %                         2*pi             2*pi
     freq_y(ii)= abs(imag(autoval(1)))/(2*pi);
     freq_z(ii)= abs(imag(autoval(2)))/(2*pi);
     freq_t(ii)= abs(imag(autoval(3)))/(2*pi);
     
     autov(:,ii) = autoval;
     lambda_old  = autoval;

     waitbar(ii/steps) 
end
close(hh)


%Plotting
figure(4)
hold on; grid on
plot(U,freq_y,'--','LineWidth',2)
plot(U,freq_z,'--','LineWidth',2)
plot(U,freq_t,'--','LineWidth',2)
xlabel('Velocity U [m/s]')
ylabel('Frequency [Hz]')
title('Frequency')
%legend('f_y','f_z','f_{\theta}')

figure(5)
hold on; grid on
plot(U,smo_y,'--','LineWidth',2)
plot(U,smo_z,'--','LineWidth',2)
plot(U,smo_t,'--','LineWidth',2)
plot(U,zeros(length(U),1),'k','LineWidth',2)
xlabel('Velocity U [m/s]')
ylabel('Damping \zeta [-]')
title('Damping')
%legend('{\zeta}_y','{\zeta}_z','{\zeta}_{\theta}')
%
figure(6)
hold on
for jj = 1:6
plot3(real(autov(jj,:)),imag(autov(jj,:)),U,'-o','LineWidth',2)
end
view(2)   %Sets the default 2-D view
grid on
xlabel('Real')
ylabel('Imag')
title('Eigenvalues')

%% TURBULENT WIND
% calcolo della U turbolenta e W turbolenta
[height,U_norm,I_u,I_w,xLu,xLw] = readvars('profiles.xlsx');

%Possibile utilizzare più cifre significative (usate 3)
%Calculus made by hand
 m_u=(U_norm(19)-U_norm(18))/(height(19)-height(18)); %slope value coef
 q_u=U_norm(19)-m_u*height(19); %origin value
 y_k_u = @(x) m_u*x + q_u;
%y_k_u = @(x) 1.72*1e-3*x + 0.87;
U_norm_deck = y_k_u(64);


 m_Iu=(I_u(19)-I_u(18))/(height(19)-height(18)); %slope value coef
 q_Iu=I_u(19)-m_Iu*height(19); %origin value
 y_k_Iu = @(x) m_Iu*x + q_Iu;
%y_k_I = @(x) -2*1e-4*x + 0.1269;
 I_u_deck = y_k_Iu(64);


 m_xL=(xLu(19)-xLu(18))/(height(19)-height(18)); %slope value coef
 q_xL=xLu(19)-m_xL*height(19); %origin value
 y_k_xL = @(x) m_xL*x + q_xL;
%y_k_xL = @(x) 1.239*x + 102.33;
 xLu_deck = y_k_xL(64);



U_ref = 45; %reference wind velocity 
T_u = xLu_deck./(U_norm_deck*U_ref);
f_sample = 10;
df = 1/600; 
f = [df:df:f_sample/2];%Nyquist sampling theory
fnu = (f.*xLu_deck)/(U_norm_deck*U_ref);
SuuN = (4*fnu)./(1+70.8*(fnu.^2)).^(5/6); % Von Karman non-dimensional power spectra
figure(7)
loglog(fnu,SuuN,'b'),grid on
xlabel('f*L/U')
ylabel('f*Suu/\sigma_u^2')
title('Von Karman non-dimensional power spectrum (u component)')

sigma_u = (I_u_deck.*(U_norm_deck*U_ref)).^2;
Suu = (SuuN.*sigma_u)./f; %spettro dimensionale
Gxx = Suu.*df;
A = sqrt(2*Gxx); %A=small u 
figure(8)
stem(f,A,'b'),grid on %Higher frequencies have lower contribution
xlabel('f [Hz]')      
ylabel('A [m/s]')
title('Power spectrum (u component)')

dt = 1/f_sample; %intervallo di campionamento
t = [0:dt:600-dt]; %time vector (10 minutes [typical length of osservation])
na = length(f); %numero armoniche
fase = rand(1,na)*2*pi; %fase random delle armoniche
Ut = ones(size(t))*U_ref; %base velocità media
for ii = 1:na
    Ut = Ut + A(ii)*sin(2*pi*f(ii)*t + fase(ii));
end
%time history plot
figure(9),
plot(t, Ut,'b'),grid on
xlabel('t [s]')
ylabel('U(t) [m/s]')
title('Time history plot of the U turbulent wind profile')
%Vettore della u piccola per un altezza di 64m ovvero H_deck
u = Ut - U_ref;

%% Faccio la stessa cosa per w

%y_k_I_w = @(x) -1.1*1e-4*x + 0.0698;
m_Iw=(I_w(19)-I_w(18))/(height(19)-height(18)); %slope value coef
q_Iw=I_w(19)-m_Iw*height(19); %origin value
y_k_Iw = @(x) m_Iw*x + q_Iw;
I_w_deck = y_k_Iw(64);
 
 
% y_k_xL_w = @(x) 0.20653*x+17.055;
m_xLw=(xLw(19)-xLw(18))/(height(19)-height(18)); %slope value coef
q_xLw=xLw(19)-m_xLw*height(19); %origin value
y_k_xLw= @(x) m_xLw*x + q_xLw;
xLw_deck = y_k_xLw(64);


U_ref = 45;
T_w = xLw_deck./(U_norm_deck*45);
f_sample = 10;
df = 1/600;
f = [df:df:f_sample/2];
fnw = f.*xLw_deck/(U_norm_deck*U_ref);
SuuN = 4*fnw./(1 + 70.8*fnw.^2).^(5/6); % Von Karman non-dimensional power spectra
figure(10)
loglog(fnw,SuuN,'b'),grid on
xlabel('f*L/U')
ylabel('f*Suu/\sigma_u^2')
title('Von Karman non-dimensional power spectrum of the w turbulent wind component')
 
sigma_w = (I_w_deck.*(U_norm_deck*U_ref)).^2;%Variance/standard deviation
Suu_w = (SuuN.*sigma_w)./f;%spettro dimensionale
Gxx_w = Suu_w.*df;
A_w = sqrt(2*Gxx_w); %A_w=small w  
figure(11)
stem(f,A_w,'b'),grid on %Higher freq.cy have lower contribution
xlabel('f [Hz]')      
ylabel('A [m/s]')
title('Dimensional power spectrum of the W turbulent wind profile')

dt = 1/f_sample; %intervallo di campionamento
t = [0:dt:600-dt]; %time vector (10 minutes [typical length of osservation])
na = length(f); %numero armoniche
fase = rand(1,na)*2*pi; %fase random delle armoniche
Wt = ones(size(t))*U_ref; %base velocità media
for ii = 1:na
    Wt = Wt + A_w(ii)*sin(2*pi*f(ii)*t + fase(ii));
end
%time history plot
figure(12)
plot(t, Wt,'b'),grid on
xlabel('t [s]')
ylabel('W(t) [m/s]')
title('Time history plot of the W turbulent wind profile')
%Vettore della u piccola per un altezza di H_deck
w = Wt - U_ref;

%% Creo la forzante di buffeting

%Non abbiamo i coefficienti aerodinamici per U = 45
%Si dispone di U(45) = 44.5570 e U(46) = 45.5696
%Per trovare il t_0 giusto per 45 m/s interpolo

m_retta = (t_statico(46) - t_statico(45))/(U(46)-U(45));
t_buff = t_statico(45) + m_retta.*(U_ref-U(45));
%Stessa cosa dell'equivalente cambio denominazione per buffeting
Cd_buff = interp1(alpha,Cd,t_buff);
Cl_buff = interp1(alpha,Cl,t_buff);
Cm_buff = interp1(alpha,Cm,t_buff);
Kd_buff = interp1(alpha,Kd,t_buff); 
Kl_buff = interp1(alpha,Kl,t_buff);
Km_buff = interp1(alpha,Km,t_buff);

f_nat = [fy; fz; ft];
om_nat = 2*pi*f_nat;
h = R./(2*M*om_nat);

%close all
for ii = [1:600] %ii=frequencies 
    R_aero = [0.5*rho*U_ref*B*2*Cd_buff   0.5*rho*U_ref*B*(Kd_buff - Cl_buff)  0;
              0.5*rho*U_ref*B*2*Cl_buff   0.5*rho*U_ref*B*(Kl_buff + Cd_buff)  0;
              0.5*rho*U_ref*B*B*2*Cm_buff 0.5*rho*U_ref*B*B*Km_buff   0.5*rho*U_ref*B*B*Km_buff*B1t];
           
    K_aero = [0 0 -0.5*rho*U_ref^2*B*Kd_buff;
              0 0 -0.5*rho*U_ref^2*B*Kl_buff;
              0 0 -0.5*rho*U_ref^2*B^2*Km_buff];
    %Total matrices   
    R = R_stru + R_aero;
    K = K_stru + K_aero;
                        
    F_bf(:,ii) = 0.5*rho*(U_ref^2)*B*[2*Cd_buff,     (Kd_buff - Cl_buff); 
                                      2*Cl_buff,     (Kl_buff + Cd_buff); 
                                      B*2*Cm_buff,        B*Km_buff]*[(A(ii)/U_ref); (A_w(ii)/U_ref)];
 
    %Vettore Omega forzante
    Omega = 2*pi*f(ii)*eye(3,3);  %la consideriamo uguale per tutti i 3 gdl 
    %F = inv(K_stru)*F_bf(:,ii);
    %a = Omega/om_nat; 
    %sqrt((1 - a.^2).^2 + (2*a*h).^2)'.*F
    
    H = -(Omega.^2)*M + 1i*Omega*R + K; %links the forces with the displacement
    x0(ii,:) = inv(H)*F_bf(:,ii); %solve equation
    %x_buff(i,:) = inv(sqrt((K - M*(Omega(:,i).^2)).^2 + (R*Omega(:,i)).^2))*F_bf;
end


figure(13)
subplot(3,1,1)
plot(f(1:240), abs(x0(1:240,1)),'r') %fino a 240 solo per avere un grafico da 0 a 0.4
hold on
plot([0.05,0.05],[0, 3],'--')
grid on
ylabel('y_buff [m]')

subplot(3,1,2)
plot(f(1:240), abs(x0(1:240,2)),'g')
hold on
plot([0.0884,0.088],[0, 0.4],'--')
grid on
ylabel('z_buff [m]')

subplot(3,1,3)
plot(f(1:240), 0.5*B*abs(x0(1:240,3)),'b')
hold on
plot([0.259,0.259],[0, 0.015],'--')
grid on
xlabel('f [Hz]')
ylabel('\theta_{buff} [m]')
%% Computing sigma displacement for all the wind 
%Displacement time history

N=600;

Pu_y= 2*x0(:,1)';
Syy=(conj(Pu_y).*Pu_y)/(4*pi*df);

Pu_z= 2*x0(:,2)';
Szz=(conj(Pu_z).*Pu_z)/(4*pi*df);

Pu_tetha=(0.5*B*2*x0(:,3))'; 
Stt=(conj(Pu_tetha).*Pu_tetha)/(4*pi*df);
% 

figure(14)
plot(f(1:N),Syy(1:N),'-*');
hold on
plot(f(1:N),Szz(1:N),'-*');
plot(f(1:N),Stt(1:N),'-*');
title('PSD')
ylabel('PSD energy [m^2/Hz]')
xlabel('Frequency [Hz]')
legend('Y-axis','Z-axis','\theta')

%% RMS

    
    hh = waitbar(0,'Computing RMS...');
    steps = 80;

for U_real=1:80

% Componnte u
    T_u = xLu_deck./(U_real);
    fn = (f.*xLu_deck)/(U_real);
    SuuN = (4*fn)./(1+70.8*(fn.^2)).^(5/6); %non-dimensional power spectral density
    sigma_u = (I_u_deck.*U_real).^2;
    Suu = (SuuN.*sigma_u)./f;
    Gxx = Suu.*df;
    A = sqrt(2*Gxx); %A=small u 
    Ut = ones(size(t))*U_real; %base velocità media
    
            for ii = 1:na
                Ut = Ut + A(ii)*sin(2*pi*f(ii)*t + fase(ii));
            end
            
    %Vettore della u piccola per un altezza di 64m ovvero H_deck
    u = Ut - U_ref;

% Componente w
    T_w = xLw_deck./(U_real);
    fn = f.*xLw_deck/(U_real);
    SuuN = 4*fn./(1 + 70.8*fn.^2).^(5/6); %non-dimensional power spectral density

    sigma_w = (I_w_deck.*(U_norm_deck*U_ref)).^2;%Variance/standard deviation
    Suu_w = (SuuN.*sigma_w)./f;%non-dimensional power spectral density
    Gxx_w = Suu_w.*df;
    A_w = sqrt(2*Gxx_w); %A_w=small w 
    Wt = ones(size(t))*U_ref; %base velocità media
    
            for ii = 1:na
                Wt = Wt + A_w(ii)*sin(2*pi*f(ii)*t + fase(ii));
            end
        
    %Vettore della u piccola per un altezza di H_deck
    w = Wt - U_ref;

    % Creo la forzante di buffeting

%Stessa cosa dell'equivalente cambio denominazione per buffeting
    Cd_buff = interp1(alpha,Cd,t_statico(U_real));
    Cl_buff = interp1(alpha,Cl,t_statico(U_real));
    Cm_buff = interp1(alpha,Cm,t_statico(U_real));
    Kd_buff = interp1(alpha,Kd,t_statico(U_real)); 
    Kl_buff = interp1(alpha,Kl,t_statico(U_real));
    Km_buff = interp1(alpha,Km,t_statico(U_real));

        for ii = 1:600
            R_aero = [0.5*rho*U_real*B*2*Cd_buff   0.5*rho*U_real*B*(Kd_buff - Cl_buff)  0;
                      0.5*rho*U_real*B*2*Cl_buff   0.5*rho*U_real*B*(Kl_buff + Cd_buff)  0;
                      0.5*rho*U_real*B*B*2*Cm_buff 0.5*rho*U_real*B*B*Km_buff   0.5*rho*U_real*B*B*Km_buff*B1t];
           
            K_aero = [0 0 -0.5*rho*U_real^2*B*Kd_buff;
              0 0 -0.5*rho*U_real^2*B*Kl_buff;
              0 0 -0.5*rho*U_real^2*B^2*Km_buff];
    %Total matrices   
            R = R_stru + R_aero;
            K = K_stru + K_aero;
                        
            F_bf(:,ii) = 0.5*rho*(U_real^2)*B*[2*Cd_buff,     (Kd_buff - Cl_buff); 
                                               2*Cl_buff,     (Kl_buff + Cd_buff); 
                                               B*2*Cm_buff,        B*Km_buff]*[(A(ii)/U_real); (A_w(ii)/U_real)];
 
    %Vettore Omega forzante
            Omega = 2*pi*f(ii)*eye(3,3);  %la consideriamo uguale per tutti i 3 gdl 
    
            H = -(Omega.^2)*M + 1i*Omega*R + K; %links the forces with the displacement
            x0(ii,:) = inv(H)*F_bf(:,ii); %solve equation
        end

% Computing the PSD for every direction

    N=600;

    Pu_y= 2*x0(:,1)';                   %Power in y direction 
    Syy=(conj(Pu_y).*Pu_y)/(4*pi*df);   %PSD in y direction

    Pu_z= 2*x0(:,2)';                   %Power in z direction
    Szz=(conj(Pu_z).*Pu_z)/(4*pi*df);   %PSD in z direction

    Pu_tetha=(0.5*B*2*x0(:,3))';                %Power in tetha direction
    Stt=(conj(Pu_tetha).*Pu_tetha)/(4*pi*df);   %PSD in tetha direction
    
%Doing average point integration method to compute sigma of each direction
%sigma formula is sigma=sqrt(Integr(PSD))in the 3 directions

    nodi = 1:1:600;
        I=0;
        h=df;  
            for i=1:((size(nodi,2))-1)    
                I=I+h.*(Syy( 1,(nodi(i+1)))+Syy (1,nodi(i)))/2;   
            end

    sigmay(U_real) = sqrt(I);
    
        I=0;
            for i=1:((size(nodi,2))-1)
                I=I+h.*(Szz( 1,(nodi(i+1)))+Szz (1,nodi(i)))/2;  
            end
            
    sigmaz(U_real) = sqrt(I);
        
        I=0;
            for i=1:((size(nodi,2))-1)
                I=I+h.*(Stt( 1,(nodi(i+1)))+Stt (1,nodi(i)))/2;
            end
            
    sigmat(U_real) = sqrt(I);

  clear I
  waitbar(U_real/steps)

end
close(hh)

    figure (15)
    
    title('RMS')
        U_real=1:80;
            plot(U_real,sigmay);
        hold on; grid on;
        xlabel('U_{average}');
        ylabel('${\sigma}$','Interpreter','latex');
            plot(U_real,sigmaz);
            plot(U_real,sigmat);
        legend('${\sigma_{y}}$','${\sigma_{z}}$','${\sigma_{t}}$','Interpreter','latex');
        title('RMS')

%% Modal Approach




