%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %------------------ ETAP 1 ------------------%
% 
% clear;
% close all;
% 
% %------------------ Ilość pomiarów ------------------%
% 
% ilosc_pomiarow = 1000;
% ilosc_a = ilosc_pomiarow/2;
% 
% %------------------ Generacja losowych wejść systemu ------------------%
% 
% u1 = randn(ilosc_pomiarow, 1);
% a1 = ones(ilosc_a, 1);
% a2 = 1.2*ones(ilosc_a, 1);
% a_zmienne = [a1; a2];    % System zmieniający się w czasie
% 
% a_stale = ones(ilosc_pomiarow, 1);     % System niezmienny w czasie
% %a_zmienne = a_stale; % Tutaj zmieniamy czy jest stałe
% 
% %------------------ Estymator jądrowy ------------------%
% 
% y1 = a_zmienne.*u1.^3;   % Definicja systemu
% %k = 0:.5:150;
% %x_wart = sin(0.1*k);
% x_wart = -10:.5:10; % Do narysowania odwrotnosci funkcji na przedziale
% pd = fitdist(x_wart', 'Kernel', 'Kernel', 'Normal');
% 
% suma_licznika = 0;
% suma_mianownika = 0;
% w = 1;
% suma_wag = 0;
% for i=length(y1):-1:1
%     K = w * pdf(pd, ((y1(i) - x_wart)/0.1));
%     uK = u1(i)* K;
%     suma_licznika = suma_licznika + uK;
%     suma_mianownika = suma_mianownika + K;
%     w = w/1.01;
%     suma_wag = suma_wag + w;
% end
% 
% miinv = suma_licznika./suma_mianownika;
% 
% %------------------ Wyznaczenie MSE ------------------%
% 
% % MSE = 0;
% % for i=1:1:length(miinv)
% %  MSE_pkt = (x_wart(i) - miinv(i))^2;
% %  MSE = MSE + MSE_pkt;
% % end
% % finalne_MSE = MSE/length(miinv);
% 
% %------------------ Generowanie spodziewanego wyjścia ------------------%
% uk = miinv;
% %y_test = 1.2.*uk.^3;
% a_k = 1.2;
% funkcja_zad = nthroot(x_wart,3)/nthroot(a_k,3);
% 
% MSE = 0;
% for i=1:1:length(miinv)
%    MSE_pkt = (funkcja_zad(i) - miinv(i))^2;
%    MSE = MSE + MSE_pkt;
% end
% finalne_MSE = MSE/length(miinv)
% 
% 
% figure(1);
% plot(x_wart,miinv);
% hold on;
% plot(x_wart,funkcja_zad);
% xlabel('yk');
% ylabel('uk');
% legend('Zbudowany model','Estymowana funkcja odwrotna');
% 
% % figure(1);
% % plot(x_wart);
% % hold on;
% % plot(y_test);
% % xlabel('k');
% % ylabel('yk');
% % legend('Oczekiwane wyjście','Uzyskane wyjście');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % ETAP 2
% 
% clear;
% close all;
% 
% % Definicja systemu
% 
% s = tf('s');
% G = 1/(s+1);
% 
% %G = c2d(G, 0.0461);
% 
% [gamma,T] = impulse(G);
% 
% % Rekurencyjny ważony estymator najmniejszych kwadratów - system niezmienny
% % w czasie
% 
% Pn = 100*eye(length(gamma));    % Definicja macierzy kowariancji
% est_gamma = zeros(length(gamma),1);   % Start algorytmu gamma0 = 0
% ff = .98;
% Phik = zeros(length(gamma),1);
% 
% i = 0;
% N_ITER = 1000;
% while(i<N_ITER)
%     if(i == 500)
%     gamma = 1.2*gamma;
%     end
%     
%     wk = randn;
%     Phik(end) = [];
%     Phik = [wk; Phik];
%     vk = sum(gamma' .* Phik');
%     zk = 0.05*randn;
%     yk = vk+zk;
%     Pn = (1/ff)*(Pn-(Pn*Phik*Phik'*Pn)/(ff+Phik'*Pn*Phik));
%     est_gamma = est_gamma+Pn*Phik*(yk-Phik'*est_gamma);
%     i = i+1;
% end
% 
% MSE = 0;
% for i=1:1:length(gamma)
%  MSE_pkt = (gamma(i) - est_gamma(i))^2;
%  MSE = MSE + MSE_pkt;
% end
% finalne_MSE = MSE/length(gamma)
% 
% x_wart = 0:1:length(gamma)-1;
% 
% figure(1);
% plot(x_wart,gamma);
% hold on;
% plot(x_wart,est_gamma);
% legend('System FIR','System wyestymowany');
% xlabel('k \epsilon [0, S]');
% ylabel('\gamma');
% 
% % Na podstawie wyznaczonego estymatora dokonujemy predykcji krok po kroku 
% 
% % k = 0:.5:150;
% % yzk = sin(0.1*k);   % Żądana funkcja 
% % 
% % Uk_pop = zeros(length(est_gamma)-1,1);
% % gamma_zero = est_gamma(1);
% % est_gamma(1) = [];
% % 
% % Wektor_sterowan = [];
% % 
% % j = 1;
% % while(j < length(k)+1)
% %     sum_gammas = sum(est_gamma' .* Uk_pop');
% %     uk = (yzk(j) - sum_gammas)/gamma_zero;
% %     Wektor_sterowan = [Wektor_sterowan; uk];
% %     Uk_pop(end) = [];
% %     Uk_pop = [uk; Uk_pop];
% %     j = j + 1;
% % end
% % 
% % pom = zeros(length(est_gamma)+1,1);
% % Y_koncowe = [];
% % 
% % c = 1;
% % while(c < length(Wektor_sterowan)+1)
% %     obecne_ster = Wektor_sterowan(c);
% %     pom(end) = [];
% %     pom = [obecne_ster; pom];
% %     vkoncowe = sum(gamma' .* pom');
% %     zkoncowe = 0.05*randn;
% %     ykoncowe = vkoncowe + zkoncowe;
% %     Y_koncowe = [Y_koncowe; ykoncowe];
% %     c = c + 1;
% % end
% % 
% % MSE = 0;
% % for i=1:1:length(Y_koncowe)
% %  MSE_pkt = (yzk(i) - Y_koncowe(i))^2;
% %  MSE = MSE + MSE_pkt;
% % end
% % finalne_MSE = MSE/length(Y_koncowe);     % Kryterium MSE do sprawdzania działania estymatora
% % 
% % figure(2);
% % plot(k, yzk);
% % hold on;
% % plot(k, Y_koncowe);
% % legend('Oczekiwane wyjście','Otrzymane wyjście');
% % xlabel('k');
% % ylabel('yk');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TESTOWA WERSJA DO BADANIA ZACHOWANIA - DLA INTUICJI 
% WERSJA TRYWIALNA - WSZYSTKIE INFORMACJE O SYSTEMIE 

% Znana postać FIR - na podstawie wyjścia wyznaczamy sygnał wewnętrzny wk
% k = 0:.5:150; % Wyjscie zadane - sinus
% y_z = sin(0.1*k);
% poprzednie_w = zeros(2,1);
% j = 1;
% tab_sygn_wew = [];
% 
% while(j < length(k)+1)
%     wkwew = y_z(j) - poprzednie_w(1) - poprzednie_w(2);
%     poprzednie_w = [wkwew;poprzednie_w];
%     poprzednie_w(end) = [];
%     tab_sygn_wew = [tab_sygn_wew;wkwew];
%     j = j + 1;
% end
% 
% i = 1;
% ak = 1;
% tab_ster = [];
% 
% while(i < length(tab_sygn_wew)+1)
%     if(i == 150)
%        ak = 1; % System stacjonarny 
%     end
%     uk = nthroot((tab_sygn_wew(i)/ak),3);
%     tab_ster = [tab_ster;uk];
%     i = i + 1;
% end
% 
% % Zadajemy sterowanie na system
% x = 1;
% Pom = zeros(3,1);
% odp_sys = [];
% x1 = -1;
% x2 = 1;
% while(x < length(tab_sygn_wew)+1)
%     if(x == 150)
%         ak = 1.2; % System stacjonarny 
%     end
%    wk = ak*tab_ster(x)^3;
%    Pom = [wk;Pom];
%    Pom(end) = [];
%    vk = sum(Pom);
%    zk = 0.05 * (x1 + ((x2-x1)*rand));
%    yk = vk + zk;
%    odp_sys = [odp_sys;yk];
%    x = x + 1;
% end
% figure(1);
% plot(odp_sys);
% hold on;
% plot(y_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) Wygenerowanie pomiarów wejścia i wyjścia po przejściu przez system

N_ITER = 10000; % ilość par pomiarowych
ak = 1; % początkowa wartość wspł. a
i = 1;
U = []; % Tablica pomiarów wejścia
Y = []; % Tablica pomiarów wyjścia
W = zeros(3,1); % Tablica pomocnicza realizująca FIR
x1 = -1; % Wejście od -1 do 1
x2 = 1;
Wektor_pomocniczy_w = [];
Szumy = [];
% Otrzymujemy dwa wektory wejść i wyjść CAŁEGO systemu

while(i < N_ITER + 1)
   if(i == 5000)
      ak = 1; % Stacjonarny/niestacjonarny 
   end
   uk = randn;
   %uk = x1 + ((x2-x1)*rand); % Generujemy pomiar z przedziału -1 do 1
   U = [U;uk]; % Dodajemy go do wektora wejść -> [starsze ...  nowsze]
   wk = ak * uk^3; % Generujemy sygnał wewnętrzny wk
   W = [wk;W]; % Dodajemy obecny pomiar do tablicy realizującej FIR
   W(end) = []; % Dawny obecny pomiar teraz staje się starszym pomiarem
   Wektor_pomocniczy_w = [Wektor_pomocniczy_w; wk];
   vk = sum(W); % Suma obecnego i dwóch starych sygnałów wk
   %zk = 0.05*randn;
   zk = 0.035/100*randn; % Zaszumienie postaci 5% sygnału wejścia
   Szumy = [Szumy;zk];
   yk = vk + zk; % Wyjście całego systemu z dodanym szumem
   Y = [Y;yk]; % Dodajemy yk do wektora wyjść -> [starsze ... nowsze]
   i = i + 1;
end

% 2) Wyznaczamy wektor wk na podstawie wyjścia yk i znajomości systemu FIR

j = 1;
poprzednie_w = zeros(2,1);
tab_sygn_wew = [];

while(j < length(Y)+1)
    wkwew = Y(j) - poprzednie_w(1) - poprzednie_w(2); % Znamy wyjście i poprzednie wejścia
    poprzednie_w = [wkwew;poprzednie_w]; % Zmiana rotacyjna 2 poprzednich wejść
    poprzednie_w(end) = [];
    tab_sygn_wew = [tab_sygn_wew; wkwew];
    j = j + 1;
end

% 3) Poszukujemy postaci sygnału wk, która zapewni nam sin na wyjściu

n = 1;
pop_w = zeros(2,1);
k = 0:.5:150;
y_z = sin(0.1*k);
tab_ster_wk = [];

while(n < length(k)+1)
   wkster = y_z(n) - pop_w(1) - pop_w(2);
   pop_w = [wkster;pop_w];
   pop_w(end) = [];
   tab_ster_wk = [tab_ster_wk;wkster];
   n = n + 1;
end

figure(1);
plot(tab_ster_wk);
xlabel('k');
ylabel('wk');

% 4) Identyfikacja na podstawie estymacji wk i oczekiwanej postaci wk

% w = 1; % Startowa postać wag
% suma_licznika = 0;
% suma_mianownika = 0;
% 
% pd = fitdist(tab_ster_wk, 'Kernel', 'Kernel', 'Normal');
% for l=length(tab_sygn_wew):-1:1
%     K = w * pdf(pd, ((tab_sygn_wew(l) - tab_ster_wk')/0.1));
%     uK = U(l)* K;
%     suma_licznika = suma_licznika + uK;
%     suma_mianownika = suma_mianownika + K;
%     w = w/1;
% end
% miinv = suma_licznika./suma_mianownika;
% 
% % a_k = 1;
% % funkcja_zad = nthroot(tab_ster_wk,3)/nthroot(a_k,3);
% % 
% % figure(1);
% % plot(tab_ster_wk,miinv);
% % hold on;
% % plot(tab_ster_wk,funkcja_zad);
% % xlabel('yk');
% % ylabel('uk');
% % legend('Zbudowany model','Estymowana funkcja odwrotna');
% 
% % % 4.5) OGRANICZENIA STEROWANIA
% % 
% % b = 1;
% % UB = 1;
% % while(b < length(miinv)+1)
% %    if(abs(miinv(b)) >= UB)
% %       if(miinv(b) < 0)
% %           miinv(b) = -UB; 
% %       end
% %       if(miinv(b) > 0)
% %          miinv(b) = UB; 
% %       end
% %    end
% %    b = b + 1; 
% % end
% % 
% %5) Zadanie sterowania na system
% 
% n = 1;
% Pom = zeros(3,1);
% odp_sys = [];
% 
% while(n < length(miinv)+1)
%     w_k = ak*miinv(n)^3;
%     Pom = [w_k;Pom];
%     Pom(end) = [];
%     v_k = sum(Pom);
%     z_k = 0.035/100*randn;
%     %z_k = 0.05*randn;
%     y_k = v_k + z_k;
%     odp_sys = [odp_sys;y_k];
%     n = n + 1;
% end
% 
% MSE = 0;
% for c=1:1:length(odp_sys)
% MSE_pkt = (y_z(c) - odp_sys(c))^2;
% MSE = MSE + MSE_pkt;
% end
% finalne_MSE = MSE/length(odp_sys) % Różnica pomiędzy estymatorem a wartością estymowaną
% Srednie_zak = sum(Szumy)/length(Szumy) % Uśrednione zakłócenie zadane w 1 etapie
% 
% figure(1);
% plot(miinv);
% xlabel('k');
% ylabel('uk');

% figure(1);
% plot(y_z);
% hold on;
% plot(odp_sys);
% xlabel('k');
% ylabel('yk');
% legend('Oczekiwane wyjście','Otrzymane wyjście')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TESTY I STARE WERSJE

% %Identyfikacja systemu na podstawie wejść wyjść i funkcji oczekiwanej
% 
% k = 0:.5:150; % Przedział na którym szukamy sin
% y_z = sin(0.1*k); % Wyjscie zadane - sinus;
% w = 1; % Startowa postać wag
% 
% suma_licznika = 0;
% suma_mianownika = 0;
% 
% pd = fitdist(y_z', 'Kernel', 'Kernel', 'Normal');
% for j=length(Y):-1:1
%     K = w * pdf(pd, ((Y(j) - y_z)/0.1));
%     uK = U(j)* K;
%     suma_licznika = suma_licznika + uK;
%     suma_mianownika = suma_mianownika + K;
%     w = w/1;
% end
% 
%  miinv = suma_licznika./suma_mianownika;
% 
% % 3) Zadajemy wyznaczone sterowanie na wejście systemu
% 
% n = 1;
% Pom = zeros(3,1);
% odp_sys = [];
% 
% while(n < length(miinv)+1) 
%     w_k = ak*miinv(n)^3;
%     Pom = [w_k;Pom];
%     Pom(end) = [];
%     v_k = sum(Pom);
%     z_k = 0.05 * (x1 + ((x2-x1)*rand));
%     y_k = v_k + z_k;
%     odp_sys = [odp_sys;y_k];
%     n = n + 1;
% end
% 
% figure(1);
% plot(y_z);
% hold on;
% plot(odp_sys);

% clear;
% close all;
% 
% % uk --> | a*u^3 | wk --> |[1,1,1]| --> vk
% 
% %%%%%%%%%%%%%%%%%
% Gamma = [1,1,1]; % Definicja podsystemu dynamicznego systemu Hammersteina
% ak = 1; % Definicja współczynnika podsystemu statycznego
% %%%%%%%%%%%%%%%%%
% N_ITER = 10000;  % Ilość par pomiarowych
% i = 0;  % Zmienna pomocnicza
% %%%%%%%%%%%%%%%%%
% U = []; % Wektor wejść systemu Hammersteina
% Y = []; % Wektor wyjść systemu Hammersteina
% W = zeros(length(Gamma),1); % Wektor pomocniczy
% %%%%%%%%%%%%%%%%%
% x1 = -1;
% x2 = 1;
% 
% % 1) Budowa wektorów pomiarowych potrzebnych w procesie identyfikacji
% 
% while(i < N_ITER)
%     if(i == N_ITER/2)
%        ak = 1.2; % Niestacjonarność badanego systemu
%     end
%     uk = x1 + ((x2-x1)*rand); % Losowy pomiar wartość oczekiwana zerowa
%     U = [uk;U];
%     wk = ak*uk^3; % Definicja podsystemu statycznego systemu Hammersteina
%     W = [wk;W]; % Przechowywanie obecnego i dwóch poprzednich pomiarów
%     W(end) = []; % Dbamy o to żeby pomiarów było równo 3
%     vk = sum(W); % Suma obecnego i dwóch poprzednich pomiarów
%     zk = 0.05 * (x1 + ((x2-x1)*rand)); % Stworzenie szumu niezależnego od procesu wejściowego
%     yk = vk + zk; % Dodanie zaszumienia do wyjścia systemu Hammersteina
%     Y = [yk;Y];
%     i = i+1;
% end
% 
% % 2) Znamy system FIR, wyznaczamy sygnał wewnętrzny wk
% 
% %%%%%%%%%%%%%%%%%
% k = 0:.5:150; % Wyjscie zadane - sinus
% y_z = sin(0.1*k);
% %%%%%%%%%%%%%%%%%
% poprzednie_w = zeros(2,1);
% j = 1;
% tab_sygn_wew = [];
% 
% while(j < length(k)+1)
%     wkwew = y_z(j) - poprzednie_w(1) - poprzednie_w(2);
%     poprzednie_w = [wkwew;poprzednie_w];
%     poprzednie_w(end) = [];
%     tab_sygn_wew = [wkwew;tab_sygn_wew];
%     j = j + 1;
% end
% 
% % 3) Znamy sygnał wewnętrzny, na podstawie pomiarów identyfikujemy blok
% % statyczny
% 
% pd = fitdist(tab_sygn_wew, 'Kernel', 'Kernel', 'Normal');
% suma_licznika = 0;
% suma_mianownika = 0;
% %%%%%%%%%%%%%%%%%
% w = 1;
% %%%%%%%%%%%%%%%%%
% 
% % Identyfikacja funkcji odwrotnej opisującej system/wyznaczenie sterowania
% for n=1:1:length(Y)
%      K = w * pdf(pd, ((Y(n) - tab_sygn_wew)/0.1));
%      uK = U(n)* K;
%      suma_licznika = suma_licznika + uK;
%      suma_mianownika = suma_mianownika + K;
%      w = w/1.1;
% end
% 
% miinv = suma_licznika./suma_mianownika; % Wyznaczona funkcja odwrotna/sterowanie
% 
% % 4) Wyznaczylismy sterowanie aby uzyskać sygnał wk, teraz podajemy je na
% % wejście systemu
% 
% d = 1;
% odp_sys = [];
% Pom = zeros(length(Gamma),1);
% % dodać niestacjonaronsoc przy zadawaniu sterowania
% while(d < length(miinv)+1)
%     w_k = ak *(miinv(d))^3; % Wyznaczenie na podstawie sterowania
%     Pom = [w_k;Pom];
%     Pom(end) = [];
%     v_k = sum(Pom);
%     z_k = 0.05 * (x1 + ((x2-x1)*rand));
%     y_k = v_k + z_k;
%     odp_sys = [y_k;odp_sys];
%     d = d+1;
% end
% 
% % Wyznaczenie MSE
% MSE = 0;
% for c=1:1:length(miinv)
% MSE_pkt = (y_z(c) - miinv(c))^2;
% MSE = MSE + MSE_pkt;
% end
% finalne_MSE = MSE/length(miinv)
% 
% % Wykresy
% 
% figure(1);
% plot(odp_sys);
% hold on;
% plot(y_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%