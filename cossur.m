%% TP communications numériques

clc
clear
close all
dbstop if error

%%
fe = 1e4;
M=4;
n_b=log2(M);
Ts=0.001;
No=1;
Ak=zeros(1,4);
Ak(1)=(1/sqrt(2))*(1+i);
Ak(2)=(1/sqrt(2))*(-1+i);
Ak(3)=(1/sqrt(2))*(1-i);
Ak(4)=(1/sqrt(2))*(-1-i);
Ns=5000;
Nb=n_b*Ns;
Nfft=512;
Te=1/fe;
Fse=Ts/Te;
bit_generes=randi([0,1],1,Nb);
Symbole=zeros(1,Ns);

g=rcosdesign(0.5,8,10);

Eg= sum(g.^2);

sigA2=1;

eb_n0_dB=0:0.5:10;
eb_n0=10.^(eb_n0_dB/10);
sigma2=sigA2*Eg./(n_b*eb_n0);

TEB=zeros(size(eb_n0));
Pb=qfunc(sqrt(2*eb_n0));

for i=1:length(eb_n0)
    error_cnt=0;
    bit_cnt=0;
    while error_cnt < 100
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
        
        ss=upsample(Symbole,Fse);
        
        
        sl=conv(ss,g);
        nl=sqrt(sigma2(i)/2)*(randn(size(sl))+1j*randn(size(sl)));
        sl=sl+nl;
        
        ga=zeros(1,length(g));
        for t=1:length(ga)
            ga(t)=conj(g(length(g)-t+1));
        end
        rl=conv(ga,sl);
        rln=rl(81:10:50160);
        symbolerecu=zeros(1,5000);
        
        for h=1:1:length(rln)
            if real(rln(h))>0 && imag(rln(h))>0
                symbolerecu(h)=0;
            elseif  real(rln(h))>0 && imag(rln(h))<0
                symbolerecu(h)=2;
            elseif  real(rln(h))<0 && imag(rln(h))>0
                symbolerecu(h)=1;
            elseif   real(rln(h))<0 && imag(rln(h))<0
                symbolerecu(h)=3;
                
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
            elseif symbolerecu(j)==3
                bitrecu(2*j-1)=1;
                bitrecu(2*j)=1;
            end
        end
        
        
        for s=1:1:10000
            
            if bit_generes(s)~= bitrecu(s)
                error_cnt=error_cnt+1;
            end
            
            
        end
        bit_cnt = bit_cnt+10000;
        
        
        
        
        
    end
    TEB(i)=error_cnt/bit_cnt;
end


semilogy(eb_n0_dB,TEB);
xlabel("Eb/No (Db)");



hold on;
grid on;
semilogy(eb_n0_dB,Pb);


legend("valeur expérimentale : TEB","valeur théorique : Pb");
hold off;

figure()
Ecart = abs(TEB-Pb);
plot(eb_n0_dB,Ecart);
grid on;
xlabel("eb n0 dB");
ylabel("ecart entre TEB et Pb");
