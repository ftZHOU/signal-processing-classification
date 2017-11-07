close all;
%clear all;

%load('analyseAvecFiltre.mat')
%load('analyseSansFiltre.mat')

figure(1);
hold on;
plot(analyseSansFiltre(2,:), '.r');
hold on;
plot(analyseAvecFiltre(2,:), '.b');
legend('sans filtre','avec filtre');
xlabel('spectre #')
ylabel('frequence maximale')