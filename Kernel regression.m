%%
clear;
clc;
close all;
%% Inicjalizacja zmiennych
N = 1000;
X_N = 4*rand(1,N)-2;
%X_N = Wzium();
Z_N = normrnd(1,1,[1,length(X_N)]);;
%Z_N = zeros(1,N);
Y_N = Nieliniowy_system_statyczny(X_N, Z_N, 1);
setGlobal(X_N,Z_N);
x = -2:0.01:2;
m = Zrob_funkcje(x,1);
jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
%jadro_box = @(x) 1/(1*sqrt(2*pi))*exp((-(x - 0).^2)/(2*1^2));
%% ZADANIE 7

gamma = 0.1;
C = @(x) 1/(pi*1*(1+((x)/1)^2));
Z_N = Generuj_liczby(C,N,-1000,1000);
Y_N = Nieliniowy_system_statyczny(X_N, Z_N, 1);
setGlobal(X_N,Z_N);
%%
mseee = MSE(m_est, m);

%% ZADANIE 1 i 2

figure(1);
plot(x,m);
hold on;
plot(X_N,Y_N,'.');
title('Zaleznosc wejscia-wyjscia systemu wraz z wygenerowana chmura pomiarow','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
%% ZADANIE 3

ZAD_3(X_N,Y_N,jadro_box,x,m);
%% ZADANIE 4

ZAD_4(X_N,Y_N,0.5,x,m); % Wykresy estymatora dla różnych jąder
%% ZADANIE 5

[h_opt,fval] = ZAD_5(); % Minimalizacja błędu
%% ZADANIE 6

m_est = Estymator(X_N, Y_N, 1, jadro_box, x);
figure(1);
plot(x,m);
hold on;
plot(x,m_est);
title('Zaleznosc wejscia-wyjscia systemu przy optymalnym wspolczynniku wygladzania i zakloceniu o niezerowej wartosci oczekiwanej','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
%% ZADANIE 6

ZAD_4(X_N,Y_N,h_opt,x,m);
%% Błąd empiryczny od h

Err = [];
for h=0.01:0.01:3
    Err(end+1) = Emp_err(10,jadro_box,h,m,x,N);
end
zakres = 0.01:0.01:3;
figure(1);
plot(zakres,Err);
title('Blad empiryczny w funkcji h','interpreter','latex');
xlabel('h','interpreter','latex');
ylabel('EmpErr(h)','interpreter','latex');
%% Błąd empiryczny od N

Err = [];
for n=50:50:5000
   Err(end+1) = Emp_err(10,jadro_box,1,m,x,n); 
end
zakres = 50:50:5000;
figure(1);
plot(zakres,Err);
title('Blad empiryczny w funkcji N przy optymalnym h','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('EmpErr(N)','interpreter','latex');
%% Szum o niezerowej wartości oczekiwanej

Z_N = normrnd(1,1,[1,N]);
Y_N = Nieliniowy_system_statyczny(X_N, Z_N, 1);
setGlobal(X_N,Z_N);

%% Cross-validation - paskudna obliczeniowo wersja

% options = optimset('Display','iter');
% [h_cros, fval_cros] = fminbnd(@CV_nie,0,5,options);
J_hN = [];

for h=0.001:0.01:2
    J_hN(end+1) = CV_nie(h);
end

zakres = 0.001:0.01:2;
figure(1);
plot(zakres,J_hN);
% hold on;
% plot(h_cros,fval_cros,'r*');
title('Przebieg wielkosci J(h) w zaleznosci od h','interpreter','latex');
xlabel('h','interpreter','latex');
ylabel('J(h)','interpreter','latex');
m_est = Estymator(X_N, Y_N, h_cros, jadro_box, x);


% plot(x,m);
% hold on;
% plot(x,m_est);
% title('Wynik identyfikacji przy optymalnym wspolczynniku h','interpreter','latex');
% xlabel('x','interpreter','latex');
% ylabel('y','interpreter','latex');

%% Kroswalidacja inteligentna wersja

% Przeszukanie zbioru danych
x0 = abs(mean(Y_N(1:10)));
mark = [];
for i=10:1:length(X_N)/10
    if(abs(mean(Y_N(i:i+10))) >= x0 * 2)
        mark(end+1) = i*10;
        x0 = abs(mean(Y_N(i:i+10)));
    elseif(abs(mean(Y_N(i:i+10))) <= x0 * 0.2)
        mark(end+1) = i*10;
        x0 = abs(mean(Y_N(i:i+10)));
    end
end



%% Funkcje i inne zabawki

function Rozklad = Wzium()
    fun = @(x) ((x>0) && (x<=1/100)) .* (50) + ((x>1/100) && (x<=1)) .* (100/198);
    U1 = rand(1,80000);
    U2 = 50*rand(1,80000);
    
    Rozklad = [];
    for i=1:1:length(U1)
       if(U2(i) <= fun(U1(i)))
          Rozklad(end+1) = U1(i); 
       end
    end
end

function err = MSE(m_est, m)
    suma = 0;
    for i=1:1:length(m)
        suma = suma + (m_est(i) - m(i))^2;
    end
    err = 1/length(m) * suma;
end

function J_h = CV(h)
    
    global X_N;
    global Z_N;
    Y_N = Nieliniowy_system_statyczny(X_N, Z_N, 1);
    n = length(X_N); % Ilość próbek
    jadro = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
    wek = [];
    
    for k=1:1:n % Leave-one-out estimator
        suma_licznik = 0;
        suma_mianownik = 0;
        for i=1:1:n
            if(i ~= k)
                K = jadro((X_N(i) - X_N(k))/h);
                suma_licznik = suma_licznik + Y_N(i)*K;
                suma_mianownik = suma_mianownik + K;
            end
        end
        wek(end+1) = (Y_N(k) - (suma_licznik/suma_mianownik))^2;
    end
    J_h = 1/n * sum(wek);
end

function J_h = CV_nie(h)
    global X_N;
    global Z_N;
    Y_N = Nieliniowy_system_statyczny(X_N, Z_N, 1);
    jadro = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
    x = -2:0.01:2;
    m_est = Estymator(X_N, Y_N, h, jadro, X_N);
    wek = [];
    for i=1:1:length(X_N)
        wek(end+1) = (Y_N(i) - m_est(i))^2;
    end
    J_h = 1/length(X_N) * sum(wek);
end

function [h_opt,fval] = ZAD_5()
    
    options = optimset('Display','iter');
    [h_opt, fval] = fminbnd(@Valid_Err,0,5,options);
    Err = [];
    for h_a=0.01:0.01:2
        Err(end+1) = Valid_Err(h_a);
    end
    zakres = 0.01:0.01:2;
    
    figure(1);
    plot(zakres,Err);
    hold on;
    plot(h_opt,fval,'r*');
    title('Blad sredniokwadratowy w funkcji h, z zaznaczonym minimum','interpreter','latex');
    xlabel('h','interpreter','latex');
    ylabel('Err','interpreter','latex');
end

function ZAD_4(X_N,Y_N,h,x,m)
    sigma = 1;
    mu = 0;
    jadro_nor = @(x) 1/(sigma*sqrt(2*pi))*exp((-(x - mu).^2)/(2*sigma^2));
    jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
    jadro_epa = @(x) ((x >= -1) && (x <= 1)) .* (3/4 *(1 - x^2)) + ((x < -1) && (x > 1)) .* 0;
    jadro_tri = @(x) ((x >= -1) && (x <= 1)) .* (70/81 * ((1 - abs(x)^3)^3))  + ((x < -1) && (x > 1)) .* 0;
    m_est1 = Estymator(X_N, Y_N, h, jadro_nor, x);
    m_est2 = Estymator(X_N, Y_N, h, jadro_box, x);
    m_est3 = Estymator(X_N, Y_N, h, jadro_epa, x);
    m_est4 = Estymator(X_N, Y_N, h, jadro_tri, x);
    
    figure(1);
    plot(x,m);
    hold on;
    plot(x,m_est1);
    hold on;
    plot(x,m_est2);
    hold on;
    plot(x,m_est3);
    hold on;
    plot(x,m_est4);
    title('Wyniki identyfikacji dla roznych funkcji jadra przy h=0.5 i a=1000','interpreter','latex');
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    legend('Oryginal','Gauss','Prost','Epa','Troj','interpreter','latex');
end

function ZAD_3(X_N,Y_N,jadro,x,m)
    h1 = 0.05;
    h2 = 0.1;
    h3 = 0.4;
    h4 = 1;
    m_est1 = Estymator(X_N, Y_N, h1, jadro, x);
    m_est2 = Estymator(X_N, Y_N, h2, jadro, x);
    m_est3 = Estymator(X_N, Y_N, h3, jadro, x);
    m_est4 = Estymator(X_N, Y_N, h4, jadro, x);
    
    figure(1);
    plot(x,m);
    hold on;
    plot(x,m_est1);
    hold on;
    plot(x,m_est2);
    hold on;
    plot(x,m_est3);
    hold on;
    plot(x,m_est4);
    title('Wyniki identyfikacji dla roznych parametrow wygladzania, jadro prostokatne','interpreter','latex');
    xlabel('x','interpreter','latex');
    ylabel('y','interpreter','latex');
    legend('Oryginal','h=0.05','h=0.1','h=0.4','h=1','interpreter','latex');
end

function Rozklad = Generuj_liczby(fun, iter, a, b)
    x = a:0.01:b;
    y = zeros(1,length(x));
    for i=1:1:length(x)
        y(i) = fun(x(i));
    end
    d = ceil(max(y));
    Rozklad = [];
    
    going = true;
    while(going)
        U1 = a + (b - a) * rand(1,1);
        U2 = d * rand(1,1);
        if(U2 <= fun(U1))
            Rozklad(end+1) = U1;
        end
        if(length(Rozklad) == iter)
            going = false;
        end
    end
end

function Err = Emp_err(L, jadro, h, y, x, N)
    sumaL = 0;
    sumaM = 0;
    M = length(x); 
    for i=1:1:L
       X_N = 4*rand(1,N)-2;
       Z_N = normrnd(1,1,[1,N]);
       Y_N = Nieliniowy_system_statyczny(X_N, Z_N, 1);
       Est = Estymator(X_N, Y_N, h, jadro, x);
       for j=1:1:M
           sumaM = sumaM + (Est(j) - y(j))^2; 
       end
       sumaL = sumaL + sumaM;
       sumaM = 0;
    end
    Err = 1/(L*M) * sumaL;
end

function setGlobal(x, z)
    global X_N;
    global Z_N;
    X_N = x;
    Z_N = z;
end

function err = Valid_Err(h)
    global X_N;
    global Z_N;
    Y_N = Nieliniowy_system_statyczny(X_N, Z_N, 1);
    jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
    x = -1:0.01:1;
    m = Zrob_funkcje(x,1);
    m_est = Estymator(X_N, Y_N, h, jadro_box, x);
    suma = 0;
    for i=1:1:length(x)
        suma = suma + (m_est(i) - m(i))^2;
    end
    err = 1/length(x) * suma;
end

function m_est = Estymator(X_N, Y_N, h, jadro, x)
    m_est = [];
    for i=1:1:length(x)
        suma_licznik = 0;
        suma_mianownik = 0;
        for j=1:1:length(X_N)
           suma_licznik = suma_licznik + (Y_N(j)*jadro((X_N(j) - x(i))/h));
           suma_mianownik = suma_mianownik + (jadro((X_N(j) - x(i))/h));
        end
        m_est(end+1) = suma_licznik/suma_mianownik;
    end
end

function Y_N = Nieliniowy_system_statyczny(X_N, Z_N, a)
    Y_N = [];
    for i=1:1:length(X_N)
        Y_N(end+1) = atan(a*X_N(i)) + Z_N(i);
    end
end

function m = Zrob_funkcje(x, a)
    m = [];
    mu = @(x) atan(a * x);
    for i=1:1:length(x)
        m(end+1) = mu(x(i));
    end
end

function Cvar = Cauchy(N)
    %cvar = @(N) tan((rand(1,N) - 0.5)*pi);
    cvar = @(N) randn(1,N)./randn(1,N);
    Cvar = cvar(N);
end
