%% etude filtre canal : 

clc 
close all
clear all

%% filtre canal à 1 trajet : 
d =2;
n=12;
nmax=(-n:n);

hn= sinc((nmax)-d);
 han= hann(length(hn));
 hn=hn.*han';

fvtool(hn);
%% filtre canal à plusieurs trajets :
d2=[1,2,3];
ao=1;
a1=1;
d1=0;
n=12;
nmax=[-n:n];
figure;

for k=1:3
hn= ao*sinc((nmax)-d1)+ a1*sinc((nmax)-d2(k));
han= hann(length(hn));
hn=hn.*han';

fvtool(hn);
figure;
plot(xcorr(hn));
title("d1 = 0 et d2 = "+k);
end


% 
% % faire autocorrelation de la réponse impulsionnelle et mesurer l'écart 
% hnbis=
% % 