clear all

s= tf('s');

% Parametry
qwe = 1;
qdiff = 0.4;
c1 = 3;
c2 = 1;
H1 = 5;
H2 = 5;
R1 = 4;
R2 = 5;
h10 = (qwe/c1)^2;
h20 = (qwe/c2)^2;
h1i = 0.2*h10;
h2i = 0.2*h20;

% Transmitancja przyrostowa
K = ((c1*H2^2)/(s*2*pi^2*R1^2*R2^2*h20^2*sqrt(h10)+pi*R2^2*h20^2*c1))/(s+(c2*H2^2)/(2*pi*R2^2*h20^2*sqrt(h20)));
[L, M] = tfdata(K, 'v');

% Stałe czasowe
% T = -1./roots(M)';
% Tmax = max(T);
T1 = (2*pi*R1^2*sqrt(h10))/c1;
T2 = (2*pi*R2^2*h20^2*sqrt(h20))/(c2*H2^2);
Tmax = max(T1, T2);

% Wzmocnienie
k = 2*sqrt(h20)/c2;

% Wykres charakterystyki statycznej h2 = f(qwe)
qgr = sqrt(H2)*c2;
qwe_static = (0:0.01:qgr)';
h2_static = (qwe_static/c2).^2;
qgr_static = h20*c2/2/sqrt(h20);
dqwe_static = (qgr_static:0.01:qgr);
dh2_static = (2*sqrt(h20)*(dqwe_static-qwe))/c2 + h20;
figure(1);
plot(qwe_static, h2_static, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
hold on;
plot(dqwe_static, dh2_static, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5);
title('Charakterystyka statyczna)');
xlabel('q');
ylabel('h');
xlim([0 qgr]);
grid on;
hold off;

% Wykres h2 oraz h2_zlinearyzowane
figure(3);
plot(h2, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
hold on;
plot(h2_zlin, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1.5);
grid on;
xlim([0 floor(20*Tmax)]);
ylim([h2i ((qwe+qdiff)/c2)^2*1.1]);
hold off;