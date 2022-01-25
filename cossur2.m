% TP communications numériques 

clc 
clear 
close all
dbstop if error
%% Variable : 

fe= 4*10^6;;
M=4;
n_b=log2(M);
Ak=zeros(1,4);
Ak(1)=(1/sqrt(2))*(1+i);
Ak(2)=(1/sqrt(2))*(-1+i);
Ak(3)=(1/sqrt(2))*(1-i);
Ak(4)=(1/sqrt(2))*(-1-i);
Ns=5000;
Nb=n_b*Ns;
Nfft=512;
Fe = 250;
Te=1/fe;
Ts=4*Te;
Fse=Ts/Te; 
bit_generes=randi([0,1],1,Nb);
Symbole=zeros(1,Ns);


%% Emetteur : 

for k=1:1:5000
    if bit_generes(2*k-1)==0 && bit_generes(2*k)==0
        Symbole(k)=Ak(1);
    elseif bit_generes(2*k-1)==0 && bit_generes(2*k)==1
        Symbole(k)=Ak(2);
    elseif bit_generes(2*k-1)==1 && bit_generes(2*k)==0
        Symbole(k)=Ak(3);
        
    elseif  bit_generes(2*k-1)==1 && bit_generes(2*k)==1
        Symbole(k)=Ak(4);
    end
    
    
end

ss=upsample(Symbole,Fse); %on echantillonne à Fse ce qui répresente environ 10 point par symboles.

g=rcosdesign(0.35,8,Fse,'sqrt');

sl=conv(ss,g);
%% canal : 
d =2;
nmax=0:25;
hn= sinc(n-12-d);
fvtool(hn);

Slavantcanal=conv(hn,sl);
var=1;
bruit= (randn(size(Slavantcanal)))+j*randn(size(Slavantcanal));
Slavantcanal=Slavantcanal+bruit;
fvtool(Slavantcanal);

ga=zeros(1,length(g));
for t=1:length(ga)-1
ga(t)=conj(g(length(g)-t));
end

rl=conv(ga,sl);
rln=rl(80:10:50160);
symbolerecu=zeros(1,5000);
%rln=downsample(rl,Fse);

for i=1:1:length(rln)-1
    if real(rln(i))>0 && imag(rln(i))>0
        symbolerecu(i)=0;
    elseif  real(rln(i))>0 && imag(rln(i))<0
        symbolerecu(i)=2;
    elseif   real(rln(i))<0 && imag(rln(i))>0
        symbolerecu(i)=1;
    elseif   real(rln(i))<0 && imag(rln(i))<0
        symbolerecu(i)=3;
    
    end
end 

bitrecu=zeros(1,10000);

for j=1:1:5000
    if symbolerecu(j)==0
        bitrecu(2*j-1)=0 ;
        bitrecu(2*j)=0;
    elseif symbolerecu(j)==1
        bitrecu(2*j-1)=0; 
        bitrecu(2*j)=1;
    elseif symbolerecu(j)==2
        bitrecu(2*j-1)=1;
        bitrecu(2*j)=0;
        
    else symbolerecu(j)==3
        bitrecu(2*j-1)=1;
        bitrecu(2*j)=1;
    end
 
  % Taux d'erreur binaire : 
  
  compteurbitfaux=0;
  nombredebits=10000;
  for i=1:1:10000
      if bit_generes(i)~= bitrecu(i)
          compteurbitfaux=compteurbitfaux+1;
      end
  end
  
  TEB = compteurbitfaux/nombredebits;
  
    
    
    
end



%% Représentation : 

abscisse=zeros(1,length(sl));

for i=1:length(sl)
    abscisse(i)=(i-1)*Te;
end
plot(abscisse,real(sl));
hold on;
grid on;
xlim([0,10*Ts-Te]);
xlabel("10Ts-Te");



abscisse2=zeros(1,length(rl));
for i=1:length(rl)
    abscisse2(i)=(i-1)*Te;
end
plot(abscisse2,real(rl));
xlim([0,10*Ts-Te]);
legend("sl(t)", "rl(t)")
title("comparaison rl(t) et sl(t) : on remarque qu'il y a un retard de 4ms")
hold off;


%periodogramme et dsp estimée :


%bennet : td pour la formule de bennet : 
Gf=fft(g,Nfft);
autocorrsl=autocorr(sl);
tfautosl=fft(autocorrsl,Nfft);
densitespectheo=(1/Ts).*tfautosl.*Gf
pas=10000/(Nfft-1);
frequence=0:pas:10000;
figure();

semilogy(frequence,abs(fftshift(densitespectheo)));
hold on;
grid on;
ylabel("densité spectrale théorique/expérimentale");
xlabel("Fréquence (Hz)");



densitespecwelch=pwelch(sl,Nfft,0,Nfft);
taille=size(densitespecwelch);

semilogy(frequence,fftshift(densitespecwelch));
legend(" densité spectrale théorique", "densité spectrale expérimentale");


hold off;

figure();
abscisse3=zeros(1,length(g));
for q=1:length(abscisse3)
    abscisse3(q)=q*Te;
end

plot(abscisse3,g);
%xlim([0.10,0.15]);
ylabel( "g cosinus sureleve");
grid on;
xlabel("Temps (s)")
title("filtre de mise en forme en cosinus surelevé");
figure();
plot(abscisse3,ga);
%xlim([0.10,0.15]);

ylabel( "ga cosinus sureleve");



%%h=h.*hann(length(h));

















 






