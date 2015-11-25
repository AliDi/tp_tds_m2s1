%TP : mod�lisation autor�gressive
%18/11/15

clear all; close all; clc

N=10; %ordre du filtre

load('bruit_structure.mat');   %chargement des donn�es du syst�me r�el
L=length(b);

autocorr_b = xcorr(b,'biased');    %calcul de l'auto-corr�lation 

%modele AR avec "Levinson-Durbin recursion"
%==========================================================================
[a_levinson, var_x]=levinson(autocorr_b(L:L+N-1),N);

%signal filtr� x avec la fonction filter et les coeff pr�c�demment trouv�s
%===========================================================================

x=sqrt(var_x)*randn(1,L); %signal � l'entr�e du filtre

y_levinson= filter(1,a_levinson,x); %signal � la sortie du filtre

autocorr_levinson=xcorr(y_levinson,'biased');

%comparaison des signaux par leurs autocorr�lation
%==========================================================================
figure
plot(temps(1:5*N),autocorr_levinson(L:L+(5*N-1)),'b');
hold on
plot(temps(1:5*N), autocorr_b(L:L+(5*N-1)),'r');
legend('R_{yy}(\tau)','R_{bb}(\tau)')
%title('Comparaison des fonction d autocorr�lation')
xlabel('\tau (s)');
ylabel('Amplitude');

%comparaison des signaux par leur contenu fr�quentiel
%===========================================================================

Te=1/fe;
f=0:fe/L:fe/2-1; %axe fr�quentiel

H_estime=0;

for n=1:N
H_estime= H_estime + a_levinson(n).*exp(-j*2*pi*n*Te*f); %fonction de transfert du filtre
end
H_estime=1./H_estime;

S_estime= (abs(H_estime)).^2.*var_x.*Te; %DSP th�orique

K=256;
[ S_estime_welch, f_welch]=pwelch(y_levinson,hanning(K),K/2,K,fe); 
[ S_reel_welch, f_welch ]= pwelch(b,hanning(K),K/2,K,fe);

figure
plot(f,10*log10(S_estime),'g')

hold on
plot(f_welch,10*log10(S_estime_welch),'b')
hold on
plot(f_welch,10*log10(S_reel_welch),'r')
ylabel('Amplitude en dB')
xlabel('Fr�quence en Hz')

legend('DSP theorique','Periodogramme de Welch sur y(t)','Periodogramme de Welch sur b(t)')

%title('Comparaison par le contenu fr�quentiel')


%recherche de l'ordre du modele le plus pertinent
%==========================================================================

for i=2:120
    [a_levinson, var_x]=levinson(autocorr_b(L:L+i-1),i);
    FPE(i)= var_x* (L+i+1)/(L-i-1);
    AIC(i) = log(var_x)+2*i/L;
end


figure
subplot(2,1,1)
plot(FPE(2:end))
ylabel('FPE')
xlabel('Ordre du mod�le')
subplot(2,1,2)
plot(AIC(2:end))
ylabel('AIC')
xlabel('Ordre du mod�le')



[a min_FPE] = min(FPE(2:end));
[a min_AIC] = min(AIC(2:end));


%recalcul de la DSP avec le bon nombre de coefficients
%==========================================================================

N=min(min_FPE,min_AIC);
[a_levinson, var_x]=levinson(autocorr_b(L:L+N),N);

H_estime=0;

for n=1:N
H_estime= H_estime + a_levinson(n).*exp(-j*2*pi*n*Te*f); %fonction de transfert du filtre
end
H_estime=1./H_estime;

S_estime= (abs(H_estime)).^2.*var_x.*Te;

figure
plot(f,10*log10(S_estime))
hold on
plot(f_welch,10*log10(S_reel_welch),'r')
title(['Comparaison des DSP avec N=' num2str(N)])
xlabel('Amplitude en dB')
ylabel('Fr�quence de Hz')

legend('DSP theorique','Periodogramme de Welch sur b(t)')






