clear ;
close all ;
clc ;

%% Initialisation des parametres
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

Eg = 0; % Energie du filtre de mise en forme ->somme des modules au carre 
for i=1:length(G)
    Eg = Eg + G(i)^2;
end

sigA2 = 1; % Variance theorique des symboles -> calcul a partir de la formule avec E(X)e

eb_n0_dB = 0:0.5:10; % Liste des Eb/N0 en dB
eb_n0 = 10.^( eb_n0_dB /10) ; % Liste des Eb/N0

sigma2 = sigA2 * Eg ./ ( n_b * eb_n0 ) ; % Variance du bruit complexe en bande de base

TEB = zeros ( size ( eb_n0 ) ); % Tableau des TEB (resultats)
Pb = qfunc ( sqrt (2* eb_n0 ) ) ; % Tableau des probabilites deerreurs theoriques = 0.5*erfc(sqrt(eb_n0))

for j = 1: length(eb_n0)
    bit_error = 0;
    bit_count = 0;
    while bit_error < 100



%% Emetteur %%

    bitemis=randi([0,3],1,Ns); %genere Ns echantillons aleatoires entre 0 et 3 (00,01,10,11)
    %1 echantillon = 2 bits           
                      
    symboleemis = pskmod(bitemis,M,pi/4,'gray'); %bit->symbole

    symboleech=upsample(symboleemis,Fse); %surechantillonnage
    
    Sl=conv2(G,symboleech); %dix echantillons = Ts en terme de temps
    %figure,plot(abs(Sl(1:50))),title("Sl avant filtre")
%% CANAL
    d=2;
    n=0:20;
    H= sinc(n-12-d).*hann(21)';
    %fvtool(H)
    H=1;
    Sl2=conv2(H,Sl);
    
    %figure,plot(abs(Sl2(1:50))),title("Sl apres filtre")
    
    nl =sqrt(sigma2(j))*(randn(size(Sl2)) + 1i*randn (size (Sl2))) ; %bruit blanc complexe
    yl = Sl2 +nl;
        
%% RECEPTEUR

    Ga = conv2(G,H); %filtre adapte
    Rg = conv2(G,Ga); %Autocorrelation entre le filtre G et le filtre adapate Ga
    
    Rh =conv2(Rg,H);
    
    retard = 0;
    max = Rh(1);
    for i=2:length(Rh)     %calcul du retard lie aux filtres
        if (Rh(i) > max)
            retard = i;
            max = Rh(i);
        end
    end
    
    rl = conv2(Ga, yl);
    
    rln = rl(retard:Fse:length(rl)-Fse); %sous-echantillonnage

    %decision 

    bitdemod = pskdemod(rln,4,pi/4,'gray'); % Symbole -> bit

    bitinit = zeros(1, n_b*length(bitemis)); % 0,1,2,3 -> 00,01,10,11 pour evaluer bit e bit
    for i=1:1:length(bitemis)
        if (bitemis(i) == 0)
            bitinit(2*i-1) = 0; bitinit(2*i) = 0;
        elseif (bitemis(i) == 1)
            bitinit(2*i-1) = 0; bitinit(2*i) = 1;
        elseif (bitemis(i) == 2)
            bitinit(2*i-1) = 1; bitinit(2*i) = 0;
        else
            bitinit(2*i-1) = 1; bitinit(2*i) = 1;
        end
    end
    
    bitfinal = zeros(1, n_b*length(bitdemod)); % 0,1,2,3 -> 00,01,10,11 pour evaluer bit e bit
    for i=1:1:length(bitdemod)
        if (bitdemod(i) == 0)
            bitfinal(2*i-1) = 0; bitfinal(2*i) = 0;
        elseif (bitdemod(i) == 1)
            bitfinal(2*i-1) = 0; bitfinal(2*i) = 1;
        elseif (bitdemod(i) == 2)
            bitfinal(2*i-1) = 1; bitfinal(2*i) = 0;
        else
            bitfinal(2*i-1) = 1; bitfinal(2*i) = 1;
        end
    end
    bit_count =0;
    bit_error =0;
    for i =1:Ns*n_b
        if(bitinit(i) ~= bitfinal(i))
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
title("evolution du TEB en fonction du SNR");


% Q2)

freq = linspace(-fech/2,fech/2, Nfft);

DSP_exp = transpose(pwelch(Sl,ones(1,Nfft),0, Nfft));

Tf_sl = fft(Sl,Nfft);
DSP_th = abs((Tf_sl)).^2;

figure(2);
semilogy(freq,fftshift(DSP_exp),'b'); %facteur 5000 d'ecart
hold on;

semilogy(freq,fftshift(real(DSP_th)),'r');
xlabel("frequence en Hz");
ylabel("DSP (s_l(t))");
title("DSP Welch en bleu DSP th en rouge");

