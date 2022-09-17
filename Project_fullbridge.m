%% Ultima parte del progetto VERSIONE 2: comparare la stabilità dell'intero ponte con 
% con quella della sezione
% Usare i file esempio , modi e le figure 

    clc
    clear all
    close all
    
    %Constants
    rho = 1.25; % air density [kg/m3]   
    B   = 35.6; % deck chord [m]
    L   = 10;    % [m] --> we have 243 nodes the bridge is 2430 m long
    
    %ERRORE NEI DATI IN INPUT CI SONO DUE SEZIONI PIù LUNGHE DELLE ALTRE
    %740 810 e quelle negative
    
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
    
    M_stru = diag([my,mz,J]);
    R_stru = diag([ry,rz,rt]);
    K_stru = diag([ky,kz,kt]);

    load('modi.mat');
%    open ('FIGURE_MODI.fig');
%    open('FEM_Izmit.fig');
 
%% Creazione delle matrici modali e delle phi per ogni grado di libertà
           
            for i=1:14   
                M_m(i,i)=parmod.m(i); % Matrice di massa modale
                K_m(i,i)=parmod.k(i); % Matrice di rigidezza modale
                R_m(i,i) = parmod.h*2*M_m(i,i)*(2*pi*parmod.freq(i));
                % Matrice di smorzamento modale formula del modello
                % sezionale generalizzata al modello completo
            end
            
    % Creazione delle 14 matrici di phi, una per ogni modo di vibrare
            
     [PHI_1]=PHI_gen(fisez,1);
     [PHI_2]=PHI_gen(fisez,2);
     [PHI_3]=PHI_gen(fisez,3);
     [PHI_4]=PHI_gen(fisez,4);
     [PHI_5]=PHI_gen(fisez,5);
     [PHI_6]=PHI_gen(fisez,6);
     [PHI_7]=PHI_gen(fisez,7);
     [PHI_8]=PHI_gen(fisez,8);
     [PHI_9]=PHI_gen(fisez,9);
     [PHI_10]=PHI_gen(fisez,10);
     [PHI_11]=PHI_gen(fisez,11);
     [PHI_12]=PHI_gen(fisez,12);
     [PHI_13]=PHI_gen(fisez,13);
     [PHI_14]=PHI_gen(fisez,14);
            
%%  Probelma statico
%Per risolvere il problema statico dobbiamo risolvere un sistema di 14
%equazioni non lineari

 % Aero Coefficents
load('IBB_polari'); %alpha is in [deg] !
Cd = DLM(:,1);    %Tutte le righe della prima colonna del file DLM: ...
                  %Contiene i valori di Cd, Cl, Cm per diversi alpha 
Cl = DLM(:,2); 
Cm = DLM(:,3);
Kd = KDLM(:,1);   %File KDLM: contiene i valori delle pendenze ...
                  %Di Cd, Cl, Cm per diversi alpha 
Kl = KDLM(:,2);
Km = KDLM(:,3);

%Generatori di phi
            
        for i=1:243
            for j=1:14
                phiy(i,j)=fisez(1,j,i);
                phiz(i,j)=fisez(2,j,i);
                phit(i,j)=fisez(3,j,i);
            end
        end

% Faccio il conto per una velocità sola

    sol_old=zeros(14,1);
    
    hh = waitbar(0,'Computing static displacemnts (full bridge)...'); 
    steps = 80;
    
        for  U=1:80
            
             options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
             qn0(:,U)=fsolve(@(qn0) NONLIN(phiy,phiz,phit,qn0,Cd,Cl,Cm,rho,B,L,alpha,K_m,U),sol_old,options);   
             sol_old=qn0(:,U);
             waitbar(U/steps) 
             
        end
      
    close(hh)
    
    
%%  Ritorno alle coordinate y z tetha 
    
    
    y_statico=(phiy)*qn0;
    z_statico=(phiz)*qn0;
    t_statico=(phit)*qn0;%Questo è in radianti
    t_statico=t_statico*(180/pi);
    
    nodi=1:243;

        figure
           
           
            subplot(1,3,1);
           
            plot(nodi,y_statico(:,4),'r');
            hold on; grid on;
            plot(nodi,y_statico(:,40),'r');
            plot(nodi,y_statico(:,80),'r-');
            title('Y-Displacement');
            xlabel('Nodes');
            ylabel('Y-axis Displacement ')

            
            subplot(1,3,2);
            plot(nodi,z_statico(:,4),'g');
            hold on; grid on;
            plot(nodi,z_statico(:,40),'g');
            plot(nodi,z_statico(:,80),'g-');
            title('Z-Displacement');
            xlabel('Nodes');
            ylabel('Z-axis Displacement')
            
            subplot(1,3,3);              
            plot(nodi,t_statico(:,4),'b');
            hold on; grid on;
            plot(nodi,y_statico(:,40),'b');
            plot(nodi,y_statico(:,80),'b-');
            title('\theta-Displacement');
            xlabel('Nodes');
            ylabel('\theta Displacement');
            hold off

        
%% Forzanti autoeccitate

    B1y = 0*B; 
    B1z = 0*B;
    B1t = 0.5*B;

        
        K_aeroM=zeros(14,14);                          
        R_aeroM=zeros(14,14);
    
       hh = waitbar(0,'Computing self exited forces (full bridge)...'); 
       steps = 80;
       
        for u=1:80
                for i=1:243 
 
                    t_0 = t_statico(:,u);       %Estrapolo le condizioni statiche per ogni velocità
                                            %theta0 -->  i Cl/k/m0;  le pend.ze Kd/l/m0
                    Cd_eq(i) = interp1(alpha,Cd,t_0(i)); % 1 riga 243 colonne
                    Cl_eq(i) = interp1(alpha,Cl,t_0(i));
                    Cm_eq(i) = interp1(alpha,Cm,t_0(i));
                    Kd_eq(i) = interp1(alpha,Kd,t_0(i));
                    Kl_eq(i) = interp1(alpha,Kl,t_0(i));
                    Km_eq(i) = interp1(alpha,Km,t_0(i));
        
                    K_aero = [0 0 -0.5*rho*u^2*B*Kd_eq(i);
                              0 0 -0.5*rho*u^2*B*Kl_eq(i);
                              0 0 -0.5*rho*u^2*B^2*Km_eq(i)];
              
                    R_aero = [0.5*rho*u*B*2*Cd_eq(i)  0.5*rho*u*B*(Kd_eq(i)-Cl_eq(i))  0;
                              0.5*rho*u*B*2*Cl_eq(i)  0.5*rho*u*B*(Kl_eq(i)+Cd_eq(i))  0;
                              0.5*rho*u*B*B*2*Cm_eq(i) 0.5*rho*u*B*B*Km_eq(i) 0.5*rho*u*B*B*Km_eq(i)*B1t];
              
                        for j=1:14
                                for ii=1:3                      
                                    phi(ii,j)=fisez(ii,j,i);                      
                                end
                        end
             
                            K_aeroM=K_aeroM+(phi')*K_aero*(phi); % phi è una 3 per 14
                            R_aeroM=R_aeroM+(phi')*R_aero*(phi); % K_aero è una 3x3
                end

                %Total matrices
                M = M_m;  %Le matrici 'stru' sono per unità di lungh.za!     
                R = R_m + R_aeroM;  %Vanno moltiplicate per B
                K = K_m + K_aeroM;
     
                % polyeig: eigenvalue/eigenvector solution
                [autovet,autoval] = polyeig(K,R,M); %Svolge ([M]lambda^2+([R]+R_aero)lambda...
                                         %                  +([K]+K_aero)phi = 0                                       
                                          
                     %Riordino gli autovalori confrontando ognuno con quello al passo
                     %precedente (in termini di modulo)
                     if u == 1
                        f_old = parmod.freq;
                        smo_old = parmod.h*ones(1,14);
                        lambda_old(1:14) = -smo_old.* (2*pi*f_old) + 1j* (2*pi*f_old) .* sqrt(1- smo_old.^2);
                        lambda_old(15:28) = -smo_old.* (2*pi*f_old) - 1j* (2*pi*f_old) .* sqrt(1- smo_old.^2);
                        %lambda_old è l'autovalore ottenuto dal sistema smorzato: 
                        %lambda_old = - ç*omega_n +- j*omega_n*radice(1 - ç^2) 
                        %           = - alpha +- j*omega_d
                     end
                     indici = 1:28;
                     for jj = 1:28
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

                     %damped frequency:    |Im(lambda)|  2*pi*f*sqrt(1-ç^2)
                     %                f_d = ---------- = ------------------
                     %                         2*pi             2*pi
                     for ii=1:14
                            smo(ii,u)= -real(autoval(ii))/norm(autoval(ii));
                            freq(ii,u)= abs(imag(autoval(ii)))/(2*pi);   % Matrice frequenze alla fine ha
                                                                            % 14 righe e 80 colonne                                                            
                     end   
                     
                autov(:,u) = autoval;
                lambda_old  = autoval;   
                
           waitbar(u/steps)    
           
        end
        
    close(hh)


%% Plotto il grafico per la stabilità
U=linspace(1,80,80);

  figure
    for i=1:14
            plot(U,smo(i,:),'--','DisplayName',['Modo',num2str(i)]);
            hold on
    end
        
        legend
    
    open ('FIGURE_MODI.fig');
    title('Full bridge stability')
    ylabel('Damping')
    xlabel('Velocita [m/s]')

 %% Plotto il grafico per la stabilità in frequency

  figure
    for i=1:14
            plot(U, freq(i,:),'DisplayName',['Frequency',num2str(i)]);
            hold on
    end
        
        legend
    
    title('Full bridge stability')
    ylabel('Frequency [Hz]')
    xlabel('Velocita [m/s]')
            
            