%TP : modélisation autorégressive
%18/11/15

clear all; close all; 

N=10; %ordre du filtre

load('bruit_structure.mat');
L=length(b);

autocorr_b = xcorr(b,'biased');

%modele AR avec "Levinson-Durbin recursion"
%==========================================================================
[a_levinson, var_x]=levinson(autocorr_b(L:L+N-1),N);

%signal filtré x avec la fonction filter et les coeff précédemment trouvés
%===========================================================================

x=sqrt(var_x)*randn(1,L); %signal à l'entrée du filtre

y_levinson= filter(1,a_levinson,x); %signal à la sortie du filtre

autocorr_levinson=xcorr(y_levinson,'biased');

%comparaison des signaux par leurs autocorrélation
%==========================================================================
figure
plot(temps(1:5*N),autocorr_levinson(L:L+(5*N-1)),'b');
hold on
plot(temps(1:5*N), autocorr_b(L:L+(5*N-1)),'r');
legend('R_{yy}(\tau)','R_{bb}(\tau)')
%title('Comparaison des fonction d autocorrélation')
xlabel('\tau (s)');
ylabel('Amplitude');

%comparaison des signaux par leur contenu fréquentiel
%===========================================================================

Te=1/fe;
f=0:fe/L:fe/2-1; %axe fréquentiel

H_estime=0;

for n=1:N
H_estime= H_estime + a_levinson(n).*exp(-j*2*pi*n*Te*f); %fonction de transfert du filtre
end
H_estime=1./H_estime;

S_estime= (abs(H_estime)).^2.*var_x.*Te; %DSP théorique

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
xlabel('Fréquence en Hz')

legend('DSP theorique','Periodogramme de Welch sur y(t)','Periodogramme de Welch sur b(t)')

%title('Comparaison par le contenu fréquentiel')


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
xlabel('Ordre du modèle')
subplot(2,1,2)
plot(AIC(2:end))
ylabel('AIC')
xlabel('Ordre du modèle')



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
ylabel('Fréquence de Hz')

legend('DSP theorique','Periodogramme de Welch sur b(t)')






