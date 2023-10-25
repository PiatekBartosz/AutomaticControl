clear all

s= tf('s');

% Parametry WALEC
R1 = 1.12838;
H1 = 5;
c1 = 4.3589;
% Parametry OSTROSŁUP
a2 = 3.59153;
b2 = 7.79612;
H2 = 19;
c2 = 3.08221;
% Parametry wejściowe
qwe = 0.5*H2;

h10 = (qwe/c1)^2;
h20 = (qwe/c2)^2;
h1i = 0;
h2i = 0.1*h20;

qdiff = sqrt(h20*1.1)*c2-qwe;

%% Stałe czasowe
T1 = 2*pi*sqrt(h10)*R1^2/c1;
T2 = 2*a2*b2*h20^(5/2)/(H2^2*c2);
T3 = 0.3*min(T1, T2);
Tmax = max(T1, T2);

% Wzmocnienie
k = 2*sqrt(h20)/c2;

%% Wzmocnienie regulatora
% kr = 2.4; % - M_max
% kr = 3.9449; % - Zapas fazy
kr = 5.2928; % - Zapas amplitudy

%% Transmitancja przyrostowa
K = k/((s*T1+1)*(s*T2+1));
[L, M] = tfdata(K, 'v');

%% Transmitancja układu otwartego
Ko = kr/(s*T3+1)*K;
[ZA, ZF] = margin(Ko);

%% Kryterium M_max
figure(1);
bodeplot((Ko/(1+Ko)), 'b');
grid on
hold on

%% Kryterium zapasu fazy
figure(2);
h1 = nyquistplot(Ko);
setoptions(h1, 'ShowFullContour', 'off', 'MagUnits', 'abs', 'PhaseUnits', 'deg');

%% Wykres uchybu w stanie ustalonym
figure(3);
plot(e, 'g', 'LineWidth', 1.2);
title('Wykres uchybu w stanie ustalonym');
xlabel('t, s');
ylabel('e(t)');
grid on;

%% Zera i bieguny układu zamkniętego
G = Ko/(1+Ko);
figure(4);
pzmap(G, 'g');
hold on;

%% Wykres h2 oraz h2_zlinearyzowane
figure(5);
plot(h2, 'b', 'LineWidth', 1.2);
grid on;
hold on;
plot(h2_zlin, 'b', 'LineWidth', 1.2);
xlim([0 floor(20*Tmax)]);
ylim([h2i ((qwe)/c2)^2*1.1]);
xlabel('t, s');
ylabel('h2');
title('Odpowiedz skokowa h2 w układzie nieliniowym oraz zlinearyzowanym');
hold off;