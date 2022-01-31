%% Initialisation
clear;
close all;
clc;

Ds = 1e6;                                       % D�bit symbole 
Ts = 1/Ds;                                      % Temps symbole
Fe=4e6;                                         % Fr�quence d'�chantillonnage 
Te=1/Fe;                                        % Temps d'�chantillonnage
Fse = Ts/Te;                                    % Facteur de sur-�chantillonnage

N= 5000;                                        % Nb de symboles par paquet
nb = 2;                                         % Nb de bits/symbole (QPSK ici)
NbBits = nb*N;                                  % Nb de bits � Tx

sigA2 = 1;                                      % Variance th�orique des symboles 
eb_n0_dB = 0:0.5:7;                            % Liste des Eb/N0 en dB
eb_n0 = 10.^( eb_n0_dB /10) ;                   % Liste des Eb/N0
G=rcosdesign(0.35,4,Fse,'sqrt');                % Filtre de mise en forme
fvtool(G)
Eg = sum(G.^2);                                 % Energie du filtre de mise en forme
sigma2 = sigA2 * Eg ./ ( nb * eb_n0 ) ;         % Variance du bruit complexe en bande de base
TEB = zeros ( size ( eb_n0 ) );                 % Tableau des TEB (r�sultats)
Pb = qfunc ( sqrt (2* eb_n0 ) ) ;               % Tableau des probabilit�s d�erreurs th�oriques = 0.5*erfc(sqrt(eb_n0))

%% Emetteur
%for j = 1: length(eb_n0)
    bit_error = 0;
    bit_count = 0;
    
    %while bit_error < 100
        sb = randi([0,1],1,NbBits);             %G�n�ration du flux binaire

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

        ssup=upsample(ss,Fse);                  %Sur�chantillonnage

        sl = conv2(G,ssup);

        %% Canal
        d=-2;
        n=0:20;
        H= sinc(n-12-d).*hann(21)';             % Canal mod�lis� par un sinus cardinal
        %H=1;
        
        fvtool(H)
        sl2=conv2(H,sl);
        
        nl =(randn(size(sl2)) + 1i*randn (size (sl2))) ; % Bruit blanc complexe

        yl= sl2 +nl;                               % Signal re�u

        %% R�cepteur
        
        R1 = conv2(G,H);                            
        Ga = conj(flip(R1));                        % Filtre adapt� � G*H
        rl = conv2(Ga, yl);
        R2 = conv2(R1,Ga);
        [~,retard] = max(R2);                      %Calcul du retard caus� par les filtres G, H et Ga

        ss_detect = rl(retard:Fse:length(rl)-Fse); %Sous-echantillonnage

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
%         bit_count =0;
%         bit_error =0;
%         for i =1:N*nb
%             if(sb(i) ~= sb_est(i))
%                 bit_error = bit_error + 1;
%             end
%             bit_count = bit_count + 1;
%         end
%     end
    %TEB(j)= bit_error/bit_count;
    BER = mean(abs(sb-sb_est));

%% Figures
figure(1)
plot(real(ss),imag(ss),'*r')
hold on;
plot(real(ss_detect),imag(ss_detect),'ob')
legend('Symboles Tx','Symboles d�tect�s')
xlabel('Partie r�elle des symboles')
ylabel('Partie imaginaire des symboles')
grid on;

% figure(2);
% semilogy(eb_n0_dB,TEB,'b');
% hold on
% semilogy(eb_n0_dB,Pb,'r');
% xlabel("E_b/N_0 en dB");
% ylabel("log(TEB)");
% title("�volution du TEB en fonction du SNR");







