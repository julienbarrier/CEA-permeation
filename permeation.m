clear all
close all
clc
format long
options = optimoptions('fmincon');
options = optimoptions(options,'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);
%on règle la tolérance pour le fit assez basse, pour que matlab trouve une
%solution correcte.

filename = './perm2/12/171211a-PET175.txt';
L = .175; %épaisseur totale du film en millimètres

delimiter = {'\t',' '};
startRow = 120;

formatSpec = '%{dd/MM/yyyy HH:mm:ss}D%f%f%f%f%f%f%[^\n\r]';
% format du fichier : date, flotants double précision (6 colonnes), séparateur de ligne

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);

temps = dataArray{:, 2};
flux = dataArray{:, 3};
temps(isnan(temps)) = []; % suppression des valeurs éventuellement vides dans le fichier
flux(isnan(flux)) = [];


clearvars filename delimiter startRow formatSpec fileID dataArray;

temps = (temps-1)*0.001558; %on met le temps en secondes. pour l'eau, les mesures sont prises avec un autre échantillonnage, remplacer par 0.002255
flux = flux-(mean(flux(1:2))); % on soustrait le bruit de fond
%figure(12)
%plot(temps,flux)

init= [mean(flux(end-40:end));.7]; % valeur d'initialisation jinf,D

Q = cumsum(flux); % intégration (méthode des rectangles)

%clearvars flux % je dégage flux de la mémoire pour ne pas encombrer inutilement le processeur

% définition de la fonction selon laquelle on fit
fonctionfit = @(a,xdata) (a(1)*L^2/a(2))*((a(2)*xdata/L^2)-(1/6)-(2/pi^2)*(-exp((-a(2)*pi^2*xdata)/(L^2))+(1/4)*exp((-4*a(2)*pi^2*xdata)/(L^2))-(1/9)*exp((-9*a(2)*pi^2*xdata)/(L^2))+(1/16)*exp((-16*a(2)*pi^2*xdata)/(L^2))-(1/25)*exp((-25*a(2)*pi^2*xdata)/(L^2))));
% fit de Q (algorithme de Levenberg-Marquardt)
opts = optimset('Algorithm', 'levenberg-marquardt');
[x,rn,res,f,opt,l,J] = lsqcurvefit(fonctionfit,init,temps,Q,[],[],opts); %fonction, paramètres, donnée x, donnée y,borne inf, borne sup, options

% renvoi des valeurs Jinf et D
Jinf = x(1)
D = x(2)

% tracé du fit pour vérification visuelle
figure(1)
plot(temps,Q,'.');
hold on
plot(temps,fonctionfit(x,temps),'linewidth',3);
xlabel('temps (h)');
legend('Q_{exp}','Q_{fit}');
ylabel('Q');


% tracé des résidus pour vérification
%figure(2)
%plot(temps,res,'.');
%title('résidus');

% tracé de la jacobienne pour vérification
%figure(3)
%plot(J(1,:),J(2,:),'linewidth',3);
%title('Jacobienne de Q');