clear;
close all;
clc;

%%Paramètres
fe=1e4;                         %%fréquence d'échantillonage
Te=1/fe;
M=4;                            %%Nbre de symboles
n_b=log2(M);                    %%Nbre de bits/symbole
Ds=1e3;                         %%Débit symbole
Ts= 1/Ds;                       %%Temps symbole
Nfft=512 ;                          %%points pour TF
Fse=Ts/Te;                          %%Facteur de sur-échantillonage (nombre d'échantillons sur une période Ts)
Ns=5000;                            %%Nombre de symboles par paquet

%%Emetteur
Sb= randi([0,1],1,Ns*n_b);   %Génération d'une séquence de bits de manière uniforme
Ss=zeros(1,Ns*n_b/2);        %Création liste stockant les symboles
k=1;
for j=1:2:Ns*n_b             % Association bits->symboles (modulation QPSK)
    if Sb(j)==0 && Sb(j+1)==0
        Ss(k)=exp(1i*pi/4);
        k=k+1;
    end
    if Sb(j)==0 && Sb(j+1)==1
        Ss(k)=exp(1i*3*pi/4);
        k=k+1;
    end
    if Sb(j)==1 && Sb(j+1)==0
        Ss(k)=exp(1i*5*pi/4);
        k=k+1;
    end
    if Sb(j)==1 && Sb(j+1)==1
        Ss(k)=exp(1i*7*pi/4);
        k=k+1;
    end
end

Ss_up=zeros(1,Ns*n_b*5);            %Sur-échantillonne le signal Ss
for i=1:Ns*n_b/2
    Ss_up((i-1)*10+1)=Ss(i);
end

G=ones(1,Fse);                       %Filtre de mise en forme porte
G2=rcosfir(0.5,4,Fse,1,'sqrt');    %Filtre de mise en forme en racine de cosinus sur-élevé

Sl=conv(G,Ss_up);           %Sortie du filtre de mise en forme porte
Sl2=conv(G2,Ss_up);         %Sortie du filtre mise en forme cos sur-élevé
