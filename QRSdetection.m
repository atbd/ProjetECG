function[ecg,necg,lp_ecg,hp_ecg,deriv_ecg,sq_ecg,moy_ecg,Fs,R_loc,QRS_width,Q_loc,S_loc]=QRSdetection(choix)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       
%              QRS detection and Rhythm Analysis                        
%                                                                       
%         QRS Beat Detection using Pan-TompKins Algorithm               
%                                                                       
% 
% 
% OUTPUTS :
% necg: signal ecg sans pic d'alimentation
% lp_ecg: signal ecg après filtre pas bas
% hp_ecg: signal ecg après filtrage pas haut
% Remarque: filtre pas bas + filtre pas haut=filtre passe bande
% deriv_ecg: signal ecg après filtre dérivateur
% sq_ecg: signal ecg après filtre quadrateur
% moy_ecg: signal ecg après filtre àmoyenne glissante
% R_loc : vecteur contenant les positions (num echantillon) des pics R detectes et localisés
% QRS_width: vecteur contenant les largeurs des QRS detectes
% Q_loc : vecteur contenant le debut des QRS detectes (marquage debut onde Q)
% S_loc : vecteur contenant la fin des QRS detectes (marquage fin onde S)
%
% remarque : QRS_width_buff=S_loc-Q_loc

 %% Loading the ECG file  %%
if (choix == 0)
 
     dossier=pwd;
    [fnam,path_scenario]=uigetfile('*.mat','Choix d''un fichier ECG');
    cmd=['cd ' path_scenario];
    eval(cmd);
    cmd=['load ' fnam(1:end-4)];
    eval(cmd);
    cmd=['cd ' dossier];
    eval(cmd);

% La fonction QRS detection a été modifiée pour faciliter les tracés et les comparaisons.
% Ainsi, le nombre donné en argument sert à choisir immédiatement l'un des
% signaux ECG sans passer par la fenêtre de choix de fichier à ouvrir
    
elseif (choix==1)
    load('ECG_PATHO1.mat');
    
elseif (choix==2)
    load('ECG_PATHO2.mat');
    
elseif (choix==3)
    load('ECG_PATHO3.mat');
    
elseif (choix==4)
    load('ECG_PATHO4.mat');

elseif (choix==5)
   % ECG de référence d'un patient sain
    load('ecg4.mat');
  %load('ECG60_2.mat');
elseif (choix == 6)
    load('ecg3.mat');
end

ecg=ecg';

necg = filtrage_encoche(ecg,Fs,60,0.01);


%% taille du fichier
Ns = length(necg);

time = [1 : Ns]/Fs;

% figure(1);
% subplot(611)
% plot(time,necg);
% ylabel('ECG brut');

% lowpass filtering 
b1 = [1 0 0 0 0 0 -2 0 0 0 0 0 1]/32;
a1 = [1 -2 1];

% figure
% freqz(b1,a1,1024,Fs)
% title('Low pass filter (0,40)')
lpecg = filter(b1, a1, necg);
lp_ecg=lpecg;

% figure(1)
% subplot(612)
% plot(time,lpecg);
% ylabel('ECG LP');

%% les différentes équations des filtres sont données dans l'article

%%% High pass filter (Allpass-lowpass) %%%
%% Low pass
b2 = zeros(1,33);
b2(1) = 1;
b2(33) = -1;
a2 = [1 -1];
% figure
% freqz(b2,a2,1024,Fs)
% title('Low pass filter ')

h2 = filter(b2,a2,lpecg);  %% lowpass

%% All pass
b3 = zeros(1,17);
b3(17) = 1;
a3 = 1;
% figure
% freqz(b3,a3,1024,Fs)
% title('All pass filter ')

h3 = filter(b3,a3,lpecg);  %% allpass

%% highpass = allpass-lowpass
p2 = h3 - h2/32;   
hp_ecg=p2;%après filtrage pass-haut

% figure(1);
% subplot(613)
% plot(time,p2);
% ylabel('ECG HP');


%% Dérivation du signal filtré
%b4 = [1 0 -1];
%b4 = b4/4;    
b4 = [1 2 0 -2 -1];
a4 = 1;

h4=filter(b4,a4,p2);    % dérivée
deriv_ecg=h4;%après dérivation

% figure(1);
% subplot(614)
% plot(time,h4);
% ylabel('Derivee');


%title('ECG after bandpass filtering and differentiation');

%%% Squaring %%%
h4 = h4.^2;
sq_ecg=h4;%après carré

% figure(1);
% subplot(615)
% plot(time,h4);
% axis tight;
% ylabel('Squared');

%% Fenetre glissante (moyenne) sur 30 echantillons
%% integration par une fenetre glissante
b5=ones(1,30)/30;
a5= 1;
op = filter(b5,a5,h4);
moy_ecg=op;

% figure(1);
% subplot(616)
% plot(time,op);
% axis tight;
% ylabel('Integrated');
% xlabel('Time in seconds');

[R_loc,QRS_width,Q_loc,S_loc]=peak_features(op,b1,b2,b3,b4,b5,necg,Fs);  

% figure;
% plot(necg,'b');
% hold on;
% stem(R_loc,necg(R_loc),'r');
% stem(Q_loc,necg(Q_loc),'g');
% stem(S_loc,necg(S_loc),'m');
% hold off;