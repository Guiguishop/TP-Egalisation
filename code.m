%% Initialisation
clear;
close all;
clc;

Ds = 1e6;                                       % Débit symbole 
Ts = 1/Ds;                                      % Temps symbole
Fe=4e6;                                         % Fréquence d'échantillonnage 
Te=1/Fe;                                        % Temps d'échantillonnage
Fse = Ts/Te;                                    % Facteur de sur-échantillonnage

N= 5000;                                        % Nb de symboles par paquet
nb = 2;                                         % Nb de bits/symbole (QPSK ici)
NbBits = nb*N;                                  % Nb de bits à Tx

 d=-2.5;
 n=0:20;
 H= sinc(n-12-d).*hann(21)';                    % Canal modélisé par un sinus cardinal

sigA2 = 1;                                      % Variance théorique des symboles 
eb_n0_dB = 0:1:10;                              % Liste des Eb/N0 en dB
eb_n0 = 10.^( eb_n0_dB /10) ;                   % Liste des Eb/N0

G=rcosdesign(0.35,4,Fse,'sqrt');                % Filtre de mise en forme
%fvtool(G)
R1 = conv2(G,H);                                %Filtre de mise en forme complet G*H                        
Ga = conj(flip(R1));                            % Filtre adapté à G*H
R2 = conv2(R1,Ga);
[Eg,retard] = max(R2);

sigma2 = sigA2 * Eg ./ ( nb * eb_n0 ) ;         % Variance du bruit complexe en bande de base
TEB = zeros ( size ( eb_n0 ) );                 % Tableau des TEB (résultats)
Pb = qfunc ( sqrt (2* eb_n0 ) ) ;               % Tableau des probabilités d’erreurs théoriques = 0.5*erfc(sqrt(eb_n0))

SNR_dB = 20;                                    % Rapport Signal/Bruit au récepteur

%% Emetteur
for j = 1: length(eb_n0)
    bit_error = 0;
    bit_count = 0;
    
    while bit_error < 100
        sb = randi([0,1],1,NbBits);             %Génération du flux binaire
        ss = zeros(1,N);
        for i=1:NbBits/nb
            tmp = sb(1,(i-1)*nb+1:i*nb);
            if tmp == [0 0]
                ss(1,i) = exp(1j*pi/4);
            elseif tmp == [0 1]
                ss(1,i) = exp(1j*3*pi/4);
            elseif tmp ==[1 1]
                ss(1,i) = exp(1j*5*pi/4);
            else
                ss(1,i) = exp(1j*7*pi/4);
            end
        end

        ssup=upsample(ss,Fse);                  %Suréchantillonnage

        sl = conv2(G,ssup);
        
        %% Canal
        d=-2.5;
        n=0:20;
        H= sinc(n-12-d).*hann(21)';             % Canal modélisé par un sinus cardinal
        %H=1;
        
        %fvtool(H)
        sl2=conv2(H,sl);
        
        Py = mean(abs(sl2).^2);               % Puissance de y
        Pbruit = Py/10^(SNR_dB/10);             % Puissance du bruit qui est ici une estimation de la variance (SNR_dB = 10log10(SNR))
        nl = sqrt(sigma2(j)/2)*(randn(size(sl2)) + 1i*randn (size (sl2))) ; % Bruit blanc complexe

        yl= sl2 +nl  ;                               % Signal reçu

        %% Récepteur    
        rl = conv2(Ga, yl);

        ss_detect = rl(retard:Fse:length(rl)-Fse); %Sous-echantillonnage
        
        %Egalisation Zero Forcing
        ed = zeros(1,length(R2));
        ed(retard)=1;
        Wzf = ed'*R2*inv(R2'*R2);
        %Wzf = ed' * pinv(R2)';
        
        ss_detect_zf = conv2(ss_detect,Wzf');
        
        ss_est = zeros(1,N);
        sb_est = zeros(size(sb));
        for i=1:N
            if real(ss_detect(1,i))>0
                if imag(ss_detect(1,i))>0
                    ss_est(1,i)=exp(1j*pi/4);
                    sb_est(1,(i-1)*nb+1:i*nb) = [0 0];
                else
                     ss_est(1,i)=exp(1j*7*pi/4);
                    sb_est(1,(i-1)*nb+1:i*nb) = [1 0];
                end
            else
                if imag(ss_detect(1,i))>0
                    ss_est(1,i)=exp(1j*3*pi/4);
                    sb_est(1,(i-1)*nb+1:i*nb) = [0 1];
                else
                    ss_est(1,i)=exp(1j*5*pi/4);
                    sb_est(1,(i-1)*nb+1:i*nb) = [1 1];
                end
            end
         end
        bit_count =0;
        bit_error =0;
        for i =1:N*nb
            if(sb(i) ~= sb_est(i))
                bit_error = bit_error + 1;
            end
            bit_count = bit_count + 1;
        end
     end
    TEB(j)= bit_error/bit_count;
end

%% Figures
figure()
plot(real(ss),imag(ss),'*r')
hold on
plot(real(ss_detect),imag(ss_detect),'ob')
legend('Symboles Tx','Symboles détectés')
xlabel('Partie réelle des symboles')
ylabel('Partie imaginaire des symboles')
grid on

figure()
semilogy(eb_n0_dB,TEB,'b');
hold on
semilogy(eb_n0_dB,Pb,'r');
xlabel("E_b/N_0 en dB");
ylabel("log(TEB)");
title("evolution du TEB en fonction du SNR");

figure()
plot(real(ss),imag(ss),'*r')
hold on
plot(real(ss_detect_zf),imag(ss_detect_zf),'ob')
legend('Symboles Tx','Symboles détectés après ZF')
xlabel('Partie réelle des symboles')
ylabel('Partie imaginaire des symboles')
grid on