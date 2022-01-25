clear;
close all;
clc;

%%Param�tres
fe=4*10^6;                         %%fr�quence d'�chantillonage
Te=1/fe;
M=4;                            %%Nbre de symboles
n_b=log2(M);                    %%Nbre de bits/symbole
Ts= 4*Te;                       %%Temps symbole
Ds=1/Ts;                         %%D�bit symbole
Nfft=512 ;                          %%points pour TF
Fse=Ts/Te;                          %%Facteur de sur-�chantillonage (nombre d'�chantillons sur une p�riode Ts)
Ns=5000;                            %%Nombre de symboles par paquet

%%Emetteur
Sb= randi([0,1],1,Ns*n_b);   %G�n�ration d'une s�quence de bits de mani�re uniforme
Ss=zeros(1,Ns*n_b/2);        %Cr�ation liste stockant les symboles
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

Ss_up=zeros(1,Ns*n_b*Fse/2);            %Sur-�chantillonne le signal Ss
for i=1:Ns*n_b/2
    Ss_up((i-1)*Fse+1)=Ss(i);
end

G=rcosdesign(0.35,8,Fse,'sqrt');    %Filtre de mise en forme en racine de cosinus sur-�lev�

Sl=conv(G,Ss_up);         %Sortie du filtre mise en forme cos sur-�lev�
