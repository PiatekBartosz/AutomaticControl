clear all;

% Parametry
k = 2;
T1 = 1;
T2 = 3;
T3 = 5;

% Zapis transmitancji
s = tf('s');
K1 = k/((1+s*T1)*(1+s*T2)*(1+s*T3));
K2 = tf(k, [T1*T2*T3 T1*T2+T1*T3+T2*T3 T1+T2+T3 1]);

% Postać zero-biegunowa
K1_zero = zpk(K1);
% Wykres zer i biegunów transmitancji
figure(1);
pzmap(K1);
title('Wykres zer i biegunów transmitancji intercji 3-go rzędu')
% Przypisanie zer, biegunów i wzmocnienia
[Z,P,W] = zpkdata(K1);
% Przypisanie w formie wektora licznika i mianownika transmitancji
[L,M] = tfdata(K1, 'v');
% Zamiana transmitancji na odpowiadające równania stanu 
% (przejście niejednoznaczne)
[A,B,C,D] = tf2ss(L, M);
% Zamiana macierzy stanów na transmitancje
% (przejście jednoznaczne)
K3 = ss(A,B,C,D);

%% 
% Odpowiedź skokowa
figure(2);
step(K1, 'r', 0:0.1:40);
title('Odpowiedź skokowa układu intercja 3-go rzędu');
ylim([0 2.5])
ylabel('A, -');
xlabel('t, s');
grid on;
% Informacje o odpowiedzi skokowej
SI = stepinfo(K1);

%% 

% Charakterystyka Nyquista
figure(3);
nyq = nyquistplot(K1);

% Wyświetlanie zapasu aplitudy liniowej
setoptions(nyq, 'ShowFullContour', 'off', 'MagUnits', 'abs')

%%

% Charakterystyka Bodego
figure(4);
bode(K1);
grid on;

% Alternatywa: bodeplot - możliwość pokazania tylko amplitudy
% bode = bodeplot(K1);
% setoptions(bode, 'MagUnits', 'abs', 'phaseVisable', 'off');

%%

% Linie pierwiastkowe
figure(5);
rlocus(K1);
kgr = margin(K1);

% linia 45 deg
% plot([-1 0], [1 0] *1) -> y = -x

%  można narywować dla konkretnego k, rlocus(K1, 'k', 2); -> czarne k=2

%  RR = rlocus(K)
%  można z RR policzyć stopień oscylacyjności

%%

% Zapas amplitudy i fazy
figure(6);
kr = 0.5*kgr;
margin(kr*K1);
[ZA, ZF] = margin(K1);

% Transmitancja główna oraz uchybowa - minreal(K, 0.1) - > dokładność przy
% skracaniu dla wyznaczania transmitancji minimalnej
G1 = minreal((kr*K1)/(1+kr*K1));
E1 = minreal(1/(1+kr*K1));

% Odpowiedź skokowa transmitancji głównej i uchybowej

%%

figure(7);
step(G1, 0:0.01:70);
hold on;
step(E1, 0:0.01:70);
grid on;
hold off;

%%

% Wskaźnik nadążania oraz regulacji

% Wskaźnik nadążania M(\omega) = G
% Wskaźnik regulacji  q(\omega) = Ge
% Jak zakłócenie w ukł zam. jest tłumiona w stosunku do otw.

% Możemy w properties wykresu zmieniać opcje takie np. jak skala 

figure(8);
bodemag(G1);

hold on
bodemag(E1);
grid on;
hold off;


% do simulinka: sim('filename.mdl') -> możemy wywoływać w pętli ze
% zmieniąjącymi się parameterami

