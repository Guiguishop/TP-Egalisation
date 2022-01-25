clear;
close all;
clc;

%%Paramètres
fe=4*10^6;                         %%fréquence d'échantillonage
Te=1/fe;
M=4;                            %%Nbre de symboles
n_b=log2(M);                    %%Nbre de bits/symbole
Ts= 4*Te;                       %%Temps symbole
Ds=1/Ts;                         %%Débit symbole
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

Ss_up=zeros(1,Ns*n_b*Fse/2);            %Sur-échantillonne le signal Ss
for i=1:Ns*n_b/2
    Ss_up((i-1)*Fse+1)=Ss(i);
end

G=rcosdesign(0.35,8,Fse,'sqrt');    %Filtre de mise en forme en racine de cosinus sur-élevé

Sl=conv(G,Ss_up);         %Sortie du filtre mise en forme cos sur-élevé
% partie canal : 
d =0;
nmax=25;
hn= sinc((0:nmax)-12-d);

fvtool(hn);
%%Récepteur
Slavantcanal=conv(hn,Sl);
var=1;
bruit= (randn(size(Slavantcanal)))+j*randn(size(Slavantcanal));
Slavantcanal=Slavantcanal+bruit;
fvtool(Slavantcanal);
Ga=zeros(1,length(G));    %Filtre adapté à un cos sur-élevé
for k=1:length(G)
    Ga(k)=G(length(G)+1-k);
end

Rl=conv(Slavantcanal,Ga);           %Sortie du filtre adapté d'un cos sur-élevé

Rl_down= Rl(length(G)-1 + (1:Fse:1 + (Ns-1)*Fse));            %Sous-échantillonne le signal Rl

Sb_final=zeros(1,Ns*n_b);                           %Association Symboles->bits par méthode du proche voisin
c=1;
for n=1:Ns*n_b/2
    if real(Rl_down(n))>0 && imag(Rl_down(n))>0        %Cas exp(1i*pi/4)
        Sb_final(c)=0;
        Sb_final(c+1)=0;
        c=c+2;
    end
    if real(Rl_down(n))<0 && imag(Rl_down(n))>0        %Cas exp(1i*3*pi/4)
        Sb_final(c)=0;
        Sb_final(c+1)=1;
        c=c+2;
    end
    if real(Rl_down(n))<=0 && imag(Rl_down(n))<=0       %Cas exp(1i*5*pi/4)
        Sb_final(c)=1;
        Sb_final(c+1)=0;
        c=c+2;
    end
    if real(Rl_down(n))>=0 && imag(Rl_down(n))<=0       %Cas exp(1i*7*pi/4)
        Sb_final(c)=1;
        Sb_final(c)=1;
        c=c+2;
    end
end

nombre_erreur= 0;
for k=1:length(Sb_final)
if (Sb(k)~=Sb_final(k))
nombre_erreur++;   
end
end
taux_erreur=
