clear all

%%%%%%%%%%%% Walce %%%%%%%%%%%%%
R1 = 1.12838;
H1 = 5;
c1 = 4.3589;

%%%%%%%%%%%% Ostroslup %%%%%%%%%%%%%
a2 = 3.59153;
b2 = 7.79612;
H2 = 19;
c2 = 3.08221;

%%%%%%%%%%%% charakterystyka statyczna %%%%%%%%%%%%%
qwe = 0.5*H2;
h10 = (qwe/c1)^2;
h20 = (qwe/c2)^2;
h1i = 0.1*h20;
h2i = 0.1*h20;
qdiff = sqrt(h20*1.1)*c2-qwe;

%%
%-----------------------------cz.1-----------------------------%
k = 2*sqrt(h20)/c2;
T1 = 2*pi*sqrt(h10)*R1^2/c1;
T2 = 2*a2*b2*h20^(5/2)/(H2^2*c2);
T3 = 0.2*min(T1, T2);

s = tf('s');
Ko=k/((1+T1*s)*(1+T2*s)*(1+T3*s));

%%

% Zapas amplitudy
amp=2; %amplituda
dk=1/amp;
eqn = @(om) pi - (atan(om * T1) + atan(om * T2) + atan(om * T3));
omm = fsolve(eqn, 0);
kra = dk*(sqrt(1 + omm^2*(T1^2))*sqrt(1 + omm^2*(T2^2))*sqrt(1 + omm^2*(T3^2)))/k

% Zapas fazy 
alf=pi/3; % zapas fazy
df=-1*pi+alf;
eqn = @(om) -df - (atan(om * T1) + atan(om * T2) + atan(om * T3));
omm = fsolve(eqn, 0);
krb = sqrt(1 + omm^2*(T1^2))*sqrt(1 + omm^2*(T2^2))*sqrt(1 + omm^2*(T3^2))/k
%%
bode(Ko);
nyquist(Ko);

kr = kra;

% % Kryterium M_max
% kr1 = 2;
% Ko1 = kr1 * Ko;
% G1 = minreal(Ko1 / (1 + Ko1));
% [M_max,fpeak] = getPeakGain(G1);
% M_max
% kr1;

% Kryterium d_K -> deltaK = 1.5
kr2=3.8571;
Ko2= kr2 * Ko;
[d_K, Pm] = margin(Ko2);
d_K
kr2;
figure
nyquist(Ko2);
bode(Ko2);

%Kryterium d_fi -> deltafi = 30stopni
kr3=3.2083;
Ko3=kr3*Ko;
[x, d_fi]=margin(Ko3);
d_fi
kr3;
figure
nyquist(Ko3);

%%%%%%%%%% wyświetlanie %%%%%%%%%%%%%%
% figure(40);
% hold on;
% step(Ko1/(1+Ko1));
% step(Ko2/(1+Ko2));
% step(Ko3/(1+Ko3));
% hold off;
% 
% figure(2);
% hold on;
% nyquist(Ko1);
% nyquist(Ko2);
% nyquist(Ko3);
% hold off;
% 
% figure(3);
% hold on;
% bode(Ko1/(1+Ko1));
% bode(Ko2/(1+Ko2));
% bode(Ko3/(1+Ko3));
% hold off;


%%
%------------------------------cz.2--------------------------------%
%%%%%%%%%% regulacja P dla kr dla delta_K %%%%%%%%%%%%%%
G2 = minreal(Ko2/(Ko2+1));
figure
step(G2)
stpinf = stepinfo(G2);
eust2 = 1 - dcgain(G2);

%%% wskaźniiki czasowe %%%
przeregulowanie_procent2 = stpinf.Overshoot;
czas_narastania2 = stpinf.RiseTime;
czas_opoznienia2 = 0; %policzyć trzeba
tr2 = stepinfo(G2, 'SettlingTimeTreshold', 0.02).TransientTime;

%%% wskaźniki częstotliwościowe %%%
[zapas_amp2, zapas_fazy2] = margin(Ko2); % <--- dla układu otwatego

% figure
% bodemag(G2);
[Max_rez2, czestotliwosc_rez2] = getPeakGain(G2);

%%% wskaźniki pierwiastkowe %%%
figure
rlocus(Ko2);
pzmap(G2);
RR2 = rlocus(Ko2, 1);
eta2 = abs(max(real(RR2)));
teta2 = max(abs(imag(RR2./real(RR2))));

zapas_amp2
zapas_fazy2
Max_rez2
przeregulowanie_procent2
tr2
eust2
eta2
teta2

%%%%%%%%%% regulacja P dla kr dla delta_fi %%%%%%%%%%%%%%
G3 = minreal(Ko3/(Ko3+1));
figure
step(G3)
stpinf1 = stepinfo(G3);
eust3 = 1 - dcgain(G3);

%%% wskaźniiki czasowe %%%
przeregulowanie_procent3 = stpinf1.Overshoot;
czas_narastania3 = stpinf1.RiseTime;
czas_opoznienia3 = 0; %policzyć trzeba
tr3 = stepinfo(G3, 'SettlingTimeTreshold', 0.02).TransientTime;

%%% wskaźniki częstotliwościowe %%%
[zapas_amp3, zapas_fazy3] = margin(Ko3);    % <--- dla układu otwatego
figure
bodemag(G3);
[Max_rez3, czestotliwosc_rez3] = getPeakGain(G3);

%%% wskaźniki pierwiastkowe %%%
figure
pzmap(G3);
RR3 = rlocus(Ko3, 1);
ni3 = abs(max(real(RR3)));
teta3 = max(abs(imag(RR3./real(RR3))));

zapas_amp3
zapas_fazy3
Max_rez3
przeregulowanie_procent3
tr3
eust3
ni3
teta3

%%
%-----------------------cz.3---------------------------------%


% Transmitancja przyrostowa
K = k/((s*T1+1)*(s*T2+1));
[L, M] = tfdata(K, 'v');

figure
data = h2.Data
plot(h2)


