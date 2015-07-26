clc
close all 
clear all



% Affichage d'un menu pour s�lectionner la partie du projet que l'on
% souhaite ex�cuter
choice = menu('Choix de la section de code � ex�cuter','1. Cr�ation d''un signal synth�tique',...
'2. Analyse spectrale','3. Filtre � encoche','4. Elimination des interf�rences sur ECG',...
'5. D�tection des complexes QRS / analyse du rythme','6. Restauration des �chantillons perdus',...
'7. Classification des pathologies');


%% 1. Cr�ation d'un signal synth�tique

% -------------------
% G�n�ration du bruit
% -------------------

% Nombre de points
f0 = 60;
fe = 200;

% Dur�e du signal � simuler
tmax = 50/f0;

N = round(fe*tmax);

% Rapport signal sur bruit
SNRdb = 100;

amp_sinus = 1;
var_bruit = amp_sinus*10.^(-SNRdb/10);

b = sqrt(var_bruit).*randn(N,1);

% ---------------------
% G�n�ration d'un sinus
% ---------------------

t = linspace(0,tmax,N)';
s = amp_sinus.*sin(2*pi*f0.*t);

s_bruite = b + s;

% On s�pare trac� et calcul
if choice == 1

    figure(1)
    subplot(2,1,1);
    plot(b)
    xlabel('t')
    ylabel('signal')
    title('Bruit blanc gaussien')
    grid on

    subplot(2,1,2);
    plot(t,s_bruite)
    grid on
    xlabel('Temps (s)')
    ylabel('Signal (V)')
    title('Sinus bruit�')
    
    % Trac� de l'histogramme
    figure(2)
    histfit(b)
    title('Histogramme de la loi normale')
    grid on 
    legend('Loi th�orique','loi simul�e')
    xlabel('VAR X')

end
%% 2. Analyse spectrale / filtrage d'un bruit d'alimentation

% --------------------------------
% 1. Periodogramme du sinus bruit�
% --------------------------------

% Z�ro-padding pour am�liorer la r�solution
N_fft = 2^(nextpow2(N) + 7);                            

% Calcul du p�riodogramme
x = (1/N)*abs(fft(s_bruite,N_fft)).^2;

% Ajout d'une fen�tre de Hanning
w = window(@hanning,N);
signal_window = w.*s_bruite;

% Calcul du p�riodogramme avec fen�tre de Hanning
x_window = (1/N)*abs(fft(signal_window,N_fft)).^2;
f = linspace(0,1,N_fft)';

% --------------------------------
% 2. Corr�logramme du sinus bruit�
% --------------------------------

% Calcul de l'autocorr�lation
autocorr_biased = xcorr(s_bruite,'biased');
% Ajout du zero-padding (donc rajout de 0 au milieu du signal)
autocorr_biased = [autocorr_biased(N:2*N-1,1);zeros(N_fft-length(autocorr_biased),1);autocorr_biased(2*N-1:-1:N+1,1)];
x_corr_biased = abs(fft(autocorr_biased));

autocorr_unbiased = xcorr(s_bruite,'unbiased');
autocorr_unbiased = [autocorr_unbiased(N:2*N-1,1);zeros(N_fft-length(autocorr_unbiased),1);autocorr_unbiased(2*N-1:-1:N+1,1)];
x_corr_unbiased = abs(fft(autocorr_unbiased));

if choice == 2

    u = -(N-1):1:(N-1);
    figure(3)
    subplot(2,1,1);
    plot(u,xcorr(s_bruite,'biased'));
    title('Trac� de l''autocorr�lation biais�e')
    xlabel('tau')
    ylabel('Valeur de l''autocorr�lation')
    grid on
    
    subplot(2,1,2);
    plot(u,xcorr(s_bruite,'unbiased'));
    title('Trac� de l''autocorr�lation non-biais�e')
    xlabel('tau')
    ylabel('Valeur de l''autocorr�lation')
    grid on
    
    % Trac� des periodogrammes
    figure(4)
    
    subplot(2,1,1)
    plot(f,10*log10(x))
    title('Trac� du p�riodogramme du sinus bruit�');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|^2 en dB')
    grid on
    
    subplot(2,1,2)
    plot(f,10*log10(x_window))
    title('Trac� du p�riodogramme du sinus bruit� avec fen�tre de Hanning');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|^2 en dB')
    grid on

    % Trac� des corr�logrammes
    figure(5)
    
    subplot(2,1,1);
    plot(f,10*log10(x_corr_biased))
    title('Trac� du corr�logramme biais� du sinus bruit�');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|')
    grid on

    subplot(2,1,2);
    plot(f,10*log10(x_corr_unbiased))
    title('Trac� du corr�logramme non-biais� du sinus bruit�');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|')
    grid on
    
    figure(6)
    plot(f,10*log10(x),'b')

    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|^2 en dB')
    title('Superposition du corr�logramme biais� et du p�riodogramme')
    grid on
    
    hold on
    plot(f,10*log10(x_corr_biased),'r--') 
    legend('P�riodogramme','Corr�logramme biais�')
    
end

%% 3. Filtre � encoche

% Calcul des constantes intervenant dans l'expression du filtre
w0 = 2*pi*f0/fe;
a1 = -2*cos(w0);

cc = hsv(9);

% Ajout d'un sinus � 10 Hz pour s'assurer qu'il ne sera pas rejet� par le
% filtre
sinus_2 = sin(2*pi*10*t);
s_bruite = s_bruite + sinus_2;

if choice == 3
    
    for epsilon = 0.1:0.1:0.9
        v = floor(10*epsilon);
        
        % Calcul des coefficients du filtre
        B = [1 a1 1];
        A = [1 (1-epsilon)*a1 (1-epsilon).^2];

        % Filtrage du signal
        y = filter(B,A,s_bruite);

        % Calcul de la fonction de transfert
        [h,w] = freqz(B,A,512);
        
        % Trac�
        figure(7)
        plot((w/pi)*fe/2,20*log10(abs(h)),'color',cc(v,:))
        title('Fonction de transfert du filtre � encoche')
        xlabel('Fr�quence (Hz)')
        ylabel('Gain (dB)')
        grid on
        hold on
    end

% Trac� du signal sinusoidal filtr� pour diff�rentes valeurs de epsilon    
% Z�ro-padding
N_fft = 2^(nextpow2(N) + 7);                            

% Calcul de la TFD
spectre_sinus = (1/N)*abs(fft(s_bruite,N_fft)).^2;
spectre_sinus_filtre = (1/N)*abs(fft(y,N_fft)).^2;     
    
legend('epsilon = 0.1','epsilon = 0.2','epsilon = 0.3','epsilon = 0.4','epsilon = 0.5',...
    'epsilon = 0.6','epsilon = 0.7','epsilon = 0.8','epsilon = 0.9','Location','SouthWest');

figure(8)
hold on;
plot(s_bruite,'red')
plot(y);
grid on
xlabel('temps (s)')
ylabel('signal')
title('signal sinuso�dal bi-ton avant et apr�s filtrage')
legend('signal original','signal filtr�')
    
% Trac� des periodogrammes
figure(9)

subplot(2,1,1)
plot(f,10*log10(spectre_sinus))
title('Trac� du p�riodogramme du sinus avant filtrage');
xlabel('Fr�quence normalis�e')
ylabel('|FFT|^2 en dB')
grid on

subplot(2,1,2)
plot(f,10*log10(spectre_sinus_filtre))
title('Trac� du p�riodogramme du sinus apr�s filtrage');
xlabel('Fr�quence normalis�e')
ylabel('|FFT|^2 en dB')
grid on    


end;


%% 2.4 Elimination des interf�rences d'alimentation pour un ECG

% Chargement d'un ECG
load('ECG60_2.mat')

% Suppression de la composante continue
ecg = ecg - mean(ecg);

N = length(ecg);

% Zero padding du signal
N_fft = 2^(nextpow2(N) + 4);   

t = linspace(0,(N-1)/Fs,N);

x = (1/N)*abs(fft(ecg,N_fft)).^2;

% Ajout d'une fenetre de Hanning
w = window(@hanning,N);
signal_window = w.*ecg;

x_window = (1/N)*abs(fft(signal_window,N_fft)).^2;
f = linspace(0,1,N_fft)';

% Calcul de l'autocorrelation
[autocorr_biased_ecg, autocorr_unbiased_ecg, x_corr_biased_ecg, x_corr_unbiased_ecg] = calcul_autocorr(ecg, N, N_fft);

ecg_filtre = zeros(N,1);

for eps = 0.1:0.3:0.9
    ecg_filtre(:,floor(10*eps)) = filtrage_encoche(ecg,60,Fs,eps)';
end

fft_ecg_filtre = ((1/N)*abs(fft(ecg_filtre,N_fft)).^2);

if choice == 4

    u = -(N-1):1:(N-1);
    figure(9)
    subplot(2,1,1);
    plot(u,xcorr(ecg,'biased'));
    title('Trac� de l''autocorr�lation biais�e de l''ECG')
    xlabel('tau')
    ylabel('Valeur de l''autocorr�lation')
    grid on
    
    subplot(2,1,2);
    plot(u,xcorr(ecg,'unbiased'));
    title('Trac� de l''autocorr�lation non-biais�e de l''ECG')
    xlabel('tau')
    ylabel('Valeur de l''autocorr�lation')
    grid on
    
    % Trac� des periodogrammes
    figure(10)
    
    subplot(2,1,1)
    plot(f,10*log10(x))
    title('Trac� du p�riodogramme de l''ECG');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|^2 en dB')
    grid on
    
    subplot(2,1,2)
    plot(f,10*log10(x_window))
    title('Trac� du p�riodogramme de l''ECG avec fen�tre de Hanning');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|^2 en dB')
    grid on

    % Trac� des corr�logrammes
    figure(11)
    
    subplot(2,1,1);
    plot(f,10*log10(x_corr_biased_ecg))
    title('Trac� du corr�logramme biais� de l''ECG');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|')
    grid on

    subplot(2,1,2);
    plot(f,10*log10(x_corr_unbiased_ecg))
    title('Trac� du corr�logramme non-biais� de l''ECG');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|')
    grid on
     
    figure(12)
    plot(t,ecg_filtre)
    title('Trac� de l''ECG filtr�');
    xlabel('temps (s)')
    ylabel('signal ECG filtre')
    grid on
    
    figure(13)
    plot(f,10*log10(fft_ecg_filtre))
    title('Trac� du p�riodogramme de l''ECG filtr� pour epsilon = 0.1, 0.4 et 0.7');
    xlabel('Fr�quence normalis�e')
    ylabel('|FFT|^2 en dB') 
    grid on
     
end


%% 3. D�tection des complexes QRS / Analyse du rythme


if choice == 5
    % Trac� des diff�rentes �tapes de la reconnaissance de complexes QRS
    [ecg,necg,lp_ecg,hp_ecg,deriv_ecg,sq_ecg,moy_ecg,Fs,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(5);
    subplot(4,1,1)
    plot(hp_ecg)
    title('Signal apr�s filtrage passe-bande')
    grid on
    
    subplot(4,1,2)
    plot(deriv_ecg)
    title('Signal apr�s filtre d�rivateur')
    grid on
    
    subplot(4,1,3)
    plot(sq_ecg)
    title('Signal apr�s filtre quadrateur')
    grid on
    
    subplot(4,1,4)
    plot(moy_ecg)
    title('Signal apr�s filtre moyenneur')
    grid on
    
    figure
    plot(ecg)
    hold on
    
    % Superposition des pics Q, R et S d�tect�s au trac� original
    plot(Q_loc,ecg(Q_loc).*ones(1,length(Q_loc)),'ko')
    
    plot(R_loc,ecg(R_loc).*ones(1,length(R_loc)),'go')
    
    plot(S_loc,ecg(S_loc).*ones(1,length(S_loc)),'mo')
    legend('ecg','Q','R','S')
    grid on
    xlabel('temps')
    ylabel('signal')
    title('Reconnaissance d''un complexe QRS')
end    


%% 4. Restauration des �chantillons perdus

if choice == 6
    
% Chargement de l'ECG
% dossier = pwd;
% [fnam,path_scenario] = uigetfile('*.mat','Choix d''un fichier ECG');
% cmd = ['cd ' path_scenario];
% eval(cmd);
% cmd = ['load ' fnam(1:end-4)];
% eval(cmd);
% cmd = ['cd ' dossier];
% eval(cmd);

% Chargement d'un ECG
load('ecg4.mat');

% Soustraction de la composante continue
ecg=ecg-mean(ecg);

% Filtrage du signal
ecg_filtre = filtrage_encoche(ecg,50,Fs,0.1);
ecg_filtre = ecg_filtre(518:691);

% Calcul du max pour normaliser le signal
max_ecg = max(ecg_filtre);

% S�lection de la zone du signal � �liminer pour tester la reconstruction..

% % Sur un complexe QRS
%x_min = 40;
%x_max = 59;

% % Sur une zone relativement uniforme
x_min = 64;
x_max = 100;

Xq = (x_min:x_max)';
X = [x_min-1 ; x_max+1];
V = ecg_filtre(X);

figure('name','Superposition');
hold on;
plot(ecg_filtre);

%ecg_filtre(x_min:x_max,1) = zeros(x_max-x_min+1,1);

% Interpolation lin�aire de la zone perdue
ecg_filtre(x_min:x_max,1) = interp1(X,V,Xq,'linear');

plot(ecg_filtre,'black');

% Ajout de zero padding
N = length(ecg_filtre);
N_fft = 2^(nextpow2(N) + 4);   
fft_ecg_filtre = zeros(N_fft,1);

% Calcul des indices du vecteur de fr�quences pour tronquer la bande 1-30
% Hz en tenant compte de la sym�trie du spectre
freqs = linspace(0,1,N_fft);
[~,a] = min(abs(freqs-(1/Fs)));
[~,b] = min(abs(freqs-(30/Fs)));
[~,c] = min(abs(freqs-(220/Fs)));
[~,d] = min(abs(freqs-(249/Fs)));

zeros1 = zeros((a-1),1);
zeros2 = zeros(c-b+1,1);
zeros3 = zeros(N_fft-(d+1),1);

% Boucle pour l'algorithme de reconstruction it�ratif
for p = 1:5
	fft_ecg_filtre = fft(ecg_filtre,N_fft); 
    
    % Troncature du spectre
    fft_ecg_filtre = [zeros1 ; fft_ecg_filtre(a:b,1) ; zeros2 ; fft_ecg_filtre(c:d,1); zeros3];
    
    % Passage � la TFD inverse pour r�cup�rer le signal dans le domaine
    % temporel
    temp = real(ifft(fft_ecg_filtre,N_fft));
    ecg_filtre(x_min:x_max) = temp(x_min:x_max);

    % On peut �galement remplacer l'int�gralit� du signal par sa TFD
    % inverse, plut�t que de remplacer uniquement la zone perdue :
    
    %ecg_filtre = real(ifft(fft_ecg_filtre,N_fft));
end

    % Trac�s des r�sultats
    %plot(real(ecg_filtre)*max_ecg/max(real(ecg_filtre)),'red')
    plot(real(ecg_filtre),'red')
    grid on
    title('Reconstruction d''une p�riode d''ECG');
    xlabel('temps')
    ylabel('signal')
    legend('original','interpol�','reconstruit')
end


%% 5. Classification des pathologies

if choice == 7

    choix = menu('S�lectionner l''ECG pathologique � �tudier :','Pathologie 1',...
        'Pathologie 2', 'Pathologie 3','Pathologie 4');

    % 5.1 Pathologie 1
    if (choix == 1)
        
        % Etude de l'ECG sain ;
        % Note : la fonction QRS d�tection a �t� modifi�e afin de prendre
        % en argument directement l'ECG que l'on souhaite �tudier, ce qui
        % �vite ainsi de passer par la fen�tre de s�lection du ficher
        [~,~,~,~,~,~,~,~,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(5);

        % Calcul de la distance entre deux ondes Q (ie la p�riode)
        periode_sain = diff(Q_loc);
        
        % Calcul de la largeur des complexes QRS
        longueurs_sain = QRS_width;

        % Calcul de la moyenne des distances entre ondes Q 
        moy_periode_sain = mean(periode_sain);
        
        % Calcul de la moyenne des longueurs QRS
        moy_longueurs_sain = mean(longueurs_sain);

        % Calcul de la variance de la p�riode
        var_periode_sain = var(periode_sain);
        
        % Calcul de la variance des complexes QRS
        var_longueurs_sain = var(longueurs_sain);

        % Calcul de la variance des distances entres ondes Q, R, et S,
        % ainsi que la variance des longueurs QRS pour un patient sain,
        % afin d'appliquer le test de Neyman Pearson
        var_Q_sain = var(Q_loc-mean(Q_loc));
        var_R_sain = var(R_loc-mean(R_loc));
        var_S_sain = var(S_loc-mean(S_loc));
        var_QRS_sain = var(QRS_width-mean(QRS_width));        
        
        
        % Etude de l'ECG malade
        [~,~,lp_ecg,hp_ecg,deriv_ecg,~,~,~,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(1);

        periode_malade = diff(Q_loc);
        longueurs_malade = QRS_width;
        
        % Comptage du nombre de pics contenus dans le signal
        % En pratique, on rel�ve tous les pics sup�rieurs � 50% de la
        % valeur de cr�te de l'ECG
        seuil = mean(hp_ecg)+0.5*(max(hp_ecg)-mean(hp_ecg));
        [X,LOCS] = findpeaks(hp_ecg,'MINPEAKHEIGHT',seuil);
        taille_peaks = length(X);
        taille_Q = length(Q_loc);
        
        % Commentaires similaires au cas pr�c�dent
        moy_periode_malade = mean(periode_malade);
        moy_longueurs_malade = mean(longueurs_malade);

        var_periode_malade = var(periode_malade);
        var_longueurs_malade = var(longueurs_malade);
        
        % Calcul du test sur les observations (somme des carr�s)
        var_Q_malade = sum((Q_loc-mean(Q_loc)).^2);
        var_R_malade = sum((R_loc-mean(R_loc)).^2);
        var_S_malade = sum((S_loc-mean(S_loc)).^2);
        var_QRS_malade = sum((QRS_width-mean(QRS_width)).^2);     
        
        % Calcul des seuils
        K1 = chi2inv(1-(10.^(-3)),length(Q_loc))*var_Q_sain;
        K2 = chi2inv(1-(10.^(-3)),length(R_loc))*var_R_sain;
        K3 = chi2inv(1-(10.^(-3)),length(S_loc))*var_S_sain;
        K4 = chi2inv(1-(10.^(-3)),length(QRS_width))*var_QRS_sain; 
        
        % V�rification (rejet ou acceptation)
        disp('Test d''hypoth�se sur la p�riode des ondes Q :')
        var_Q_malade > K1
        disp('Test d''hypoth�se sur la p�riode des ondes R :')
        var_R_malade > K2
        disp('Test d''hypoth�se sur la p�riode des ondes S :')
        var_S_malade > K3
        disp('Test d''hypoth�se sur la largeur des ondes QRS :')
        var_QRS_malade > K4
        
    end
    
    % 5.2 Pathologie 2
    if (choix == 2)
        
        % Etude de l'ECG sain ;
        % Note : la fonction QRS d�tection a �t� modifi�e afin de prendre
        % en argument directement l'ECG que l'on souhaite �tudier, ce qui
        % �vite ainsi de passer par la fen�tre de s�lection du ficher
        [~,~,~,~,~,~,~,~,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(5);

        % Calcul de la distance entre deux ondes Q (ie la p�riode)
        periode_sain = diff(Q_loc);
        
        % Calcul de la largeur des complexes QRS
        longueurs_sain = QRS_width;

        % Calcul de la moyenne des distances entre ondes Q 
        moy_periode_sain = mean(periode_sain);
        
        % Calcul de la moyenne des longueurs QRS
        moy_longueurs_sain = mean(longueurs_sain);

        % Calcul de la variance de la p�riode
        var_periode_sain = var(periode_sain);
        
        % Calcul de la variance des complexes QRS
        var_longueurs_sain = var(longueurs_sain);

        % Calcul de la variance des distances entres ondes Q, R, et S,
        % ainsi que la variance des longueurs QRS pour un patient sain,
        % afin d'appliquer le test de Neyman Pearson
        var_Q_sain = var(Q_loc-mean(Q_loc));
        var_R_sain = var(R_loc-mean(R_loc));
        var_S_sain = var(S_loc-mean(S_loc));
        var_QRS_sain = var(QRS_width-mean(QRS_width));        
        
        
        % Etude de l'ECG malade
        [~,~,lp_ecg,hp_ecg,deriv_ecg,~,~,~,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(2);

        periode_malade = diff(Q_loc);
        longueurs_malade = QRS_width;
        
        % Comptage du nombre de pics contenus dans le signal
        % En pratique, on rel�ve tous les pics sup�rieurs � 50% de la
        % valeur de cr�te de l'ECG
        seuil = mean(hp_ecg)+0.5*(max(hp_ecg)-mean(hp_ecg));
        [X,LOCS] = findpeaks(hp_ecg,'MINPEAKHEIGHT',seuil);
        taille_peaks = length(X);
        taille_Q = length(Q_loc);
        
        % Commentaires similaires au cas pr�c�dent
        moy_periode_malade = mean(periode_malade);
        moy_longueurs_malade = mean(longueurs_malade);

        var_periode_malade = var(periode_malade);
        var_longueurs_malade = var(longueurs_malade);
        
        % Calcul du test sur les observations (somme des carr�s)
        var_Q_malade = sum((Q_loc-mean(Q_loc)).^2);
        var_R_malade = sum((R_loc-mean(R_loc)).^2);
        var_S_malade = sum((S_loc-mean(S_loc)).^2);
        var_QRS_malade = sum((QRS_width-mean(QRS_width)).^2);     
        
        % Calcul des seuils
        K1 = chi2inv(1-(10.^(-3)),length(Q_loc))*var_Q_sain;
        K2 = chi2inv(1-(10.^(-3)),length(R_loc))*var_R_sain;
        K3 = chi2inv(1-(10.^(-3)),length(S_loc))*var_S_sain;
        K4 = chi2inv(1-(10.^(-3)),length(QRS_width))*var_QRS_sain; 
        
        % V�rification (rejet ou acceptation)
        disp('Test d''hypoth�se sur la p�riode des ondes Q :')
        var_Q_malade > K1
        disp('Test d''hypoth�se sur la p�riode des ondes R :')
        var_R_malade > K2
        disp('Test d''hypoth�se sur la p�riode des ondes S :')
        var_S_malade > K3
        disp('Test d''hypoth�se sur la largeur des ondes QRS :')
        var_QRS_malade > K4        
        
        
    end
    
    % 5.3 Pathologie 3
    if (choix == 3)
        
        % Etude de l'ECG sain ;
        % Note : la fonction QRS d�tection a �t� modifi�e afin de prendre
        % en argument directement l'ECG que l'on souhaite �tudier, ce qui
        % �vite ainsi de passer par la fen�tre de s�lection du ficher
        [~,~,~,~,~,~,~,~,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(5);

        % Calcul de la distance entre deux ondes Q (ie la p�riode)
        periode_sain = diff(Q_loc);
        
        % Calcul de la largeur des complexes QRS
        longueurs_sain = QRS_width;

        % Calcul de la moyenne des distances entre ondes Q 
        moy_periode_sain = mean(periode_sain);
        
        % Calcul de la moyenne des longueurs QRS
        moy_longueurs_sain = mean(longueurs_sain);

        % Calcul de la variance de la p�riode
        var_periode_sain = var(periode_sain);
        
        % Calcul de la variance des complexes QRS
        var_longueurs_sain = var(longueurs_sain);

        % Calcul de la variance des distances entres ondes Q, R, et S,
        % ainsi que la variance des longueurs QRS pour un patient sain,
        % afin d'appliquer le test de Neyman Pearson
        var_Q_sain = var(Q_loc-mean(Q_loc));
        var_R_sain = var(R_loc-mean(R_loc));
        var_S_sain = var(S_loc-mean(S_loc));
        var_QRS_sain = var(QRS_width-mean(QRS_width));        
        
        
        % Etude de l'ECG malade
        [~,~,lp_ecg,hp_ecg,deriv_ecg,~,~,~,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(3);

        periode_malade = diff(Q_loc);
        longueurs_malade = QRS_width;
        
        % Comptage du nombre de pics contenus dans le signal
        % En pratique, on rel�ve tous les pics sup�rieurs � 50% de la
        % valeur de cr�te de l'ECG
        seuil = mean(hp_ecg)+0.5*(max(hp_ecg)-mean(hp_ecg));
        [X,LOCS] = findpeaks(hp_ecg,'MINPEAKHEIGHT',seuil);
        taille_peaks = length(X);
        taille_Q = length(Q_loc);
        
        % Commentaires similaires au cas pr�c�dent
        moy_periode_malade = mean(periode_malade);
        moy_longueurs_malade = mean(longueurs_malade);

        var_periode_malade = var(periode_malade);
        var_longueurs_malade = var(longueurs_malade);
        
        % Calcul du test sur les observations (somme des carr�s)
        var_Q_malade = sum((Q_loc-mean(Q_loc)).^2);
        var_R_malade = sum((R_loc-mean(R_loc)).^2);
        var_S_malade = sum((S_loc-mean(S_loc)).^2);
        var_QRS_malade = sum((QRS_width-mean(QRS_width)).^2);     
        
        % Calcul des seuils
        K1 = chi2inv(1-(10.^(-3)),length(Q_loc))*var_Q_sain;
        K2 = chi2inv(1-(10.^(-3)),length(R_loc))*var_R_sain;
        K3 = chi2inv(1-(10.^(-3)),length(S_loc))*var_S_sain;
        K4 = chi2inv(1-(10.^(-3)),length(QRS_width))*var_QRS_sain; 
        
        % V�rification (rejet ou acceptation)
        disp('Test d''hypoth�se sur la p�riode des ondes Q :')
        var_Q_malade > K1
        disp('Test d''hypoth�se sur la p�riode des ondes R :')
        var_R_malade > K2
        disp('Test d''hypoth�se sur la p�riode des ondes S :')
        var_S_malade > K3
        disp('Test d''hypoth�se sur la largeur des ondes QRS :')
        var_QRS_malade > K4      
        
    end
    
     % 5.4 Pathologie 4
     if (choix == 4)
        
        % Etude de l'ECG sain ;
        % Note : la fonction QRS d�tection a �t� modifi�e afin de prendre
        % en argument directement l'ECG que l'on souhaite �tudier, ce qui
        % �vite ainsi de passer par la fen�tre de s�lection du ficher
        [~,~,~,~,~,~,~,~,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(5);

        % Calcul de la distance entre deux ondes Q (ie la p�riode)
        periode_sain = diff(Q_loc);
        
        % Calcul de la largeur des complexes QRS
        longueurs_sain = QRS_width;

        % Calcul de la moyenne des distances entre ondes Q 
        moy_periode_sain = mean(periode_sain);
        
        % Calcul de la moyenne des longueurs QRS
        moy_longueurs_sain = mean(longueurs_sain);

        % Calcul de la variance de la p�riode
        var_periode_sain = var(periode_sain);
        
        % Calcul de la variance des complexes QRS
        var_longueurs_sain = var(longueurs_sain);

        % Calcul de la variance des distances entres ondes Q, R, et S,
        % ainsi que la variance des longueurs QRS pour un patient sain,
        % afin d'appliquer le test de Neyman Pearson
        var_Q_sain = var(Q_loc-mean(Q_loc));
        var_R_sain = var(R_loc-mean(R_loc));
        var_S_sain = var(S_loc-mean(S_loc));
        var_QRS_sain = var(QRS_width-mean(QRS_width));        
        
        
        % Etude de l'ECG malade
        [~,~,lp_ecg,hp_ecg,deriv_ecg,~,~,~,R_loc,QRS_width,Q_loc,S_loc] = QRSdetection(4);

        periode_malade = diff(Q_loc);
        longueurs_malade = QRS_width;
        
        % Comptage du nombre de pics contenus dans le signal
        % En pratique, on rel�ve tous les pics sup�rieurs � 50% de la
        % valeur de cr�te de l'ECG
        seuil = mean(hp_ecg)+0.5*(max(hp_ecg)-mean(hp_ecg));
        [X,LOCS] = findpeaks(hp_ecg,'MINPEAKHEIGHT',seuil);
        taille_peaks = length(X);
        taille_Q = length(Q_loc);
        
        % Commentaires similaires au cas pr�c�dent
        moy_periode_malade = mean(periode_malade);
        moy_longueurs_malade = mean(longueurs_malade);

        var_periode_malade = var(periode_malade);
        var_longueurs_malade = var(longueurs_malade);
        
        % Calcul du test sur les observations (somme des carr�s)
        var_Q_malade = sum((Q_loc-mean(Q_loc)).^2);
        var_R_malade = sum((R_loc-mean(R_loc)).^2);
        var_S_malade = sum((S_loc-mean(S_loc)).^2);
        var_QRS_malade = sum((QRS_width-mean(QRS_width)).^2);     
        
        % Calcul des seuils
        K1 = chi2inv(1-(10.^(-3)),length(Q_loc))*var_Q_sain;
        K2 = chi2inv(1-(10.^(-3)),length(R_loc))*var_R_sain;
        K3 = chi2inv(1-(10.^(-3)),length(S_loc))*var_S_sain;
        K4 = chi2inv(1-(10.^(-3)),length(QRS_width))*var_QRS_sain; 
        
        % V�rification (rejet ou acceptation)
        disp('Test d''hypoth�se sur la p�riode des ondes Q :')
        var_Q_malade > K1
        disp('Test d''hypoth�se sur la p�riode des ondes R :')
        var_R_malade > K2
        disp('Test d''hypoth�se sur la p�riode des ondes S :')
        var_S_malade > K3
        disp('Test d''hypoth�se sur la largeur des ondes QRS :')
        var_QRS_malade > K4         
        
        
     end
% -----------------------------------------------
% Reconnaissance de pathologie, version empirique
% -----------------------------------------------

% Apr�s �tude, on est arriv� � l'arbre de d�cision suivant :

     if (abs((moy_periode_malade-moy_periode_sain)/moy_periode_sain) < 0.3)
        msgbox('Le patient souffre de la pathologie 3','Resultat du test empirique')      
     else    
        if (abs((var_longueurs_malade-var_longueurs_sain)/var_longueurs_sain) < 0.3)
            msgbox('Le patient souffre de la pathologie 4','Resultat du test empirique')      
        else
            if (abs((taille_peaks - taille_Q)/taille_Q) < 0.3)
                msgbox('Le patient souffre de la pathologie 1','Resultat du test empirique')      
            else
                msgbox('Le patient souffre de la pathologie 2','Resultat du test empirique')     
            end    
        end 
     end



end    













































