clear all
close all
clc

load('bruit_structure.mat');  
%signal utile b
%frequence d'echantillonnage fe
% axe_temps temps

%==========================================================================
%Q1: Estimateur de la fonction d'autocorrelation
Er = xcorr(b,'biased');

%==========================================================================
%Q2: modèle autoregressif à 10 coefs
ordre_modele = 20;                      %nb de points pris pour calculer les coefs du filtre
[coefs e] = levinson(Er(length(temps):length(temps)+ordre_modele));      % calcul des coefs du filtre
bruit =  sqrt(e)*randn(1,length(b));                                            % calcul d'un bruit
modele = filter(1,coefs,bruit);                                        % filtrage du bruit

%comparaison des fct d'autocorr du modele et du signal
figure(1)
plot(-length(temps)+1:length(temps)-1, xcorr(b));
hold on
plot(-length(temps)+1:length(temps)-1, xcorr(modele),'-r');
legend('signal','modele');
title('comparaison des autocorrelations');

%comparaison des signaux temporelles
figure(2)
plot(temps,b);
hold on
plot(temps, modele,'-r');
legend('signal','modele');
title('comparaison des signaux');
xlabel('temps');

Ntfd = length(temps);
f = (0:(Ntfd-1))*(fe/Ntfd/2);
% figure(3)
% plot(f,abs(fft(b)))
% hold on
% plot(f,abs(fft(modele)),'-r')
% xlim([0 fe/2]);
% legend('signal','modele');
% xlabel('frequence en Hz');
% ylabel('fft');
% title('comparaison des fft');
% 

% le modele et le signal ne se supperposent pas dans le domaine temporelle
% cependant les autocorrelation montrent que les 2 signaux ou les mêmes
% caractéristiques statistiques: on a donc bien une modelisation du signal
% les signaux ont le mm contenu frequentiel
%==========================================================================

% calcul de la DSP theorique du modele   (diapo 43)
V =e;

    S = 0;
for n=1:ordre_modele                                  % somme sur les coefs
     S = S +coefs(n) * exp(-i*2*pi*n*f/fe); 
end    
DSP_modele = V ./ abs(S).^2;


%affichage de la DSP theorique du modele
figure(4)
plot(f, 10*log10(DSP_modele/fe));
title('DSP theorique du modele calcul diapo 43');
ylabel('10 log10(DSP)');
xlabel('frequence en Hz');

% Calcul et affichages de periodogramme de welch
K = 256;
[Sth, axef1 ] = pwelch(modele, hanning(K), K/2, K,fe); 
[Sb,  axef2 ] = pwelch(b, hanning(K), K/2, K,fe);

figure(5)
plot(axef2,10*log10(Sb));
hold on
plot(axef1,10*log10(Sth),'-r');
legend('signal','modele');
title(['Comparaison des periodogramme de welch pour K = ' num2str(K) ' points']);
ylabel('10 log10(periodogramme)');
xlabel('frequence en Hz');


%====================================================================================
N=500;
FPE = zeros(1,N);
AIC = zeros(1,N);

for x=1:1:N
    [coefs V] = aryule(b,x);
    FPE(x) = V*(Ntfd+x+1)/(Ntfd-x-1);
    AIC(x) = log(V)+2*x/Ntfd;
    V
end

figure(6)
plot(FPE);
hold on
plot(AIC,'-r');
title('Comparaison des 2 criteres');
xlabel('ordre du modele');
ylabel('erreur');
legend('FPE','AIC');

%=============================================================
[x y] = min(AIC);
ordre_modele = y;
[coefs V] = aryule(b,y);

for n=1:ordre_modele                                  % somme sur les coefs
     S = S +coefs(n) * exp(-i*2*pi*n*f/fe); 
end    
DSP_modele = V ./ abs(S).^2;

figure(10)
plot(f, 10*log10(DSP_modele/fe));
title(['DSP theorique du modele avec ordre opti o=' num2str(ordre_modele)]);
ylabel('10 log10(DSP)');
xlabel('frequence en Hz');