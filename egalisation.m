clear ;
close all ;
clc ;

%% Initialisation des param�tres
fech=4*10^6;
Te=1/fech;
Ds= 10^6;
Ts=1/Ds;
Fse=Ts/Te;

M=4;
n_b=log2(M);
Ak=[(-1-1j)/sqrt(2); (-1+1j)/sqrt(2); (1-1j)/sqrt(2); (1+1j)/sqrt(2)];

Ns=5000;
Nc=Ns*Fse;
Nfft=512;

%Choix du filtre
G=rcosdesign(0.35,4,Fse,'sqrt');  %filtre de mise en forme

Eg = 0; % Energie du filtre de mise en forme ->somme des modules au carr� 
for i=1:length(G)
    Eg = Eg + G(i)^2;
end

sigA2 = 1; % Variance th�orique des symboles -> calcul a partir de la formule avec E(X)�

eb_n0_dB = 0:0.5:10; % Liste des Eb/N0 en dB
eb_n0 = 10.^( eb_n0_dB /10) ; % Liste des Eb/N0

sigma2 = sigA2 * Eg ./ ( n_b * eb_n0 ) ; % Variance du bruit complexe en bande de base

TEB = zeros ( size ( eb_n0 ) ); % Tableau des TEB (r�sultats)
Pb = qfunc ( sqrt (2* eb_n0 ) ) ; % Tableau des probabilit�s d�erreurs th�oriques = 0.5*erfc(sqrt(eb_n0))

for j = 1: length(eb_n0)
    bit_error = 0;
    bit_count = 0;
    while bit_error < 100



%% Emetteur %%

    Sb=randi([0,3],1,Ns); %g�n�re Ns �chantillons al�atoires entre 0 et 3 (00,01,10,11)
    %1 �chantillon = 2 bits           
                      
    Ss = pskmod(Sb,M,pi/4,'gray'); %bit->symbole

    Ss2=upsample(Ss,Fse); %sur�chantillonnage
    
    Sl=conv2(G,Ss2); %dix �chantillons = Ts en terme de temps
    %figure,plot(abs(Sl(1:50))),title("Sl avant filtre")
%% CANAL
    d=2;
    n=0:20;
    H= sinc(n-12-d).*hann(21)';
    %fvtool(H)
    Sl2=conv2(H,Sl);
    
    %figure,plot(abs(Sl2(1:50))),title("Sl apr�s filtre")
    
    nl =sqrt(sigma2(j))*(randn(size(Sl2)) + 1i*randn (size (Sl2))) ; %bruit blanc complexe
    yl = Sl2 +nl;
        
%% RECEPTEUR

    Ga = conv2(G,H); %filtre adapt�
    Rg = conv2(G,Ga); %Autocorr�lation entre le filtre G et le filtre adapat� Ga
    
    Rh =conv2(Rg,H);
    
    retard = 0;
    max = Rh(1);
    for i=2:length(Rh)     %calcul du retard li� aux filtres
        if (Rh(i) > max)
            retard = i;
            max = Rh(i);
        end
    end
    
    rl = conv2(Ga, yl);
    
    rln = rl(retard:Fse:length(rl)); %sous-echantillonnage

    %d�cision 

    bn = pskdemod(rln,4,pi/4,'gray'); % Symbole -> bit

    Sb2 = zeros(1, n_b*length(Sb)); % 0,1,2,3 -> 00,01,10,11 pour �valuer bit � bit
    for i=1:1:length(Sb)
        if (Sb(i) == 0)
            Sb2(2*i-1) = 0; Sb2(2*i) = 0;
        elseif (Sb(i) == 1)
            Sb2(2*i-1) = 0; Sb2(2*i) = 1;
        elseif (Sb(i) == 2)
            Sb2(2*i-1) = 1; Sb2(2*i) = 0;
        else
            Sb2(2*i-1) = 1; Sb2(2*i) = 1;
        end
    end
    
    bn2 = zeros(1, n_b*length(bn)); % 0,1,2,3 -> 00,01,10,11 pour �valuer bit � bit
    for i=1:1:length(bn)
        if (bn(i) == 0)
            bn2(2*i-1) = 0; bn2(2*i) = 0;
        elseif (bn(i) == 1)
            bn2(2*i-1) = 0; bn2(2*i) = 1;
        elseif (bn(i) == 2)
            bn2(2*i-1) = 1; bn2(2*i) = 0;
        else
            bn2(2*i-1) = 1; bn2(2*i) = 1;
        end
    end
    bit_count =0;
    bit_error =0;
    for i =1:Ns*n_b
        if(Sb2(i) ~= bn2(i))
            bit_error = bit_error + 1;
        end
        bit_count = bit_count + 1;
    end
    end
    TEB(j)= bit_error/bit_count;
end

%%Affichage
figure();
semilogy(eb_n0_dB,TEB,'b');
hold on
semilogy(eb_n0_dB,Pb,'r');
xlabel("E_b/N_0 en dB");
ylabel("log(TEB)");
title("�volution du TEB en fonction du SNR");


% Q2)

freq = linspace(-fech/2,fech/2, Nfft);

DSP_exp = transpose(pwelch(Sl,ones(1,Nfft),0, Nfft));

Tf_sl = fft(Sl,Nfft);
DSP_th = abs((Tf_sl)).^2;

figure(2);
semilogy(freq,fftshift(DSP_exp),'b'); %facteur 5000 d'�cart
hold on;

semilogy(freq,fftshift(real(DSP_th)),'r');
xlabel("fr�quence en Hz");
ylabel("DSP (s_l(t))");
title("DSP Welch en bleu DSP th en rouge");

