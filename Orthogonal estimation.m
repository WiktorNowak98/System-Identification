%%
clear;
clc;

%%
x = -4:0.01:4;
a = 1;
N = 2500;

m = Zrob_funkcje_kotek(x,a);
%m = Zrob_funkcje_exp(x);
X_N = 2*pi*rand(1,N)-pi;
%X_N = 6*rand(1,N)-2;
Z_N = normrnd(0,1,[1,N]);
%Z_N = zeros(1,N);
Y_N = Nieliniowy_system_statyczny(X_N, Z_N, a);
%Y_N = System_statyczny_exp(X_N, Z_N);
setGlobal(X_N, Y_N);

%% Zakłócenie cauchy'ego
C = @(x) 1/(pi*1*(1+((x)/1)^2));
Z_N = Generuj_liczby(C,N,-1000,1000);
Y_N = Nieliniowy_system_statyczny(X_N, Z_N, a);
setGlobal(X_N, Y_N);

%% Charakterystyka systemu wraz z chmurą pomiarów
figure(1);
plot(x,m)
hold on;
plot(X_N,Y_N,'.');
title('Charakterystyka badanego systemu wraz z wygenerowana chmura pomiarow','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');

%% Sprawozdanie - zad 4
L1 = 2;
L2 = 3;
L3 = 4;
L4 = 5;

m1 = Estymator_ortogonalny(X_N, Y_N, x, L1);
m2 = Estymator_ortogonalny(X_N, Y_N, x, L2);
m3 = Estymator_ortogonalny(X_N, Y_N, x, L3);
m4 = Estymator_ortogonalny(X_N, Y_N, x, L4);

% m1 = Estymator_Fourier(X_N, Y_N, x, L1, 7);
% m2 = Estymator_Fourier(X_N, Y_N, x, L2, 7);
% m3 = Estymator_Fourier(X_N, Y_N, x, L3, 7);
% m4 = Estymator_Fourier(X_N, Y_N, x, L4, 7);

subplot(2,2,1);
plot(x,m);
hold on;
plot(x,m1);
title('Ortogonalny estymator funkcji regresji L=2','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
subplot(2,2,2);
plot(x,m);
hold on;
plot(x,m2);
title('Ortogonalny estymator funkcji regresji L=3','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
subplot(2,2,3);
plot(x,m);
hold on;
plot(x,m3);
title('Ortogonalny estymator funkcji regresji L=4','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
subplot(2,2,4);
plot(x,m);
hold on;
plot(x,m4);
title('Ortogonalny estymator funkcji regresji L=5','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');

%% Estymator ortogonalny
L = 100;
m_est = Estymator_ortogonalny(X_N, Y_N, x, L);

figure(1);
plot(x,m);
hold on;
plot(x,real(m_est));

%% Błąd średniokwadratowy

options = optimset('Display','iter');
[L_opt, fval] = fminbnd(@ValidError,1,45,options);

Err = [];
for L=1:1:45
    Err(end+1) = ValidError(L);
end
zakres = 1:1:45;

subplot(2,1,1);
plot(zakres,Err);
hold on;
plot(L_opt,fval,'r*');
title('Minimalizacja bledu walidacji zaleznego od wspolczynnika L','interpreter','latex');
xlabel('L','interpreter','latex');
ylabel('Err','interpreter','latex');

est_opt = Estymator_ortogonalny(X_N, Y_N, x, L_opt);
subplot(2,1,2);
plot(x,m);
hold on;
plot(x,est_opt);
title('Ortogonalny estymator funkcji regresji - zminimalizowany blad','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');

%% Baza fouriera
L = 10;
est = Estymator_Fourier(X_N, Y_N, x, L, 7);

figure(1);
plot(x,m);
hold on;
plot(x,real(est));
title('Ortogonalny estymator funkcji regresji - baza Fouriera','interpreter','latex');
xlabel('x','interpreter','latex');
ylabel('y','interpreter','latex');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Estymator                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vErr = ValidError(L)
    global X_N;
    global Y_N;
    x = linspace(-3,3,100);
    m = Zrob_funkcje_kotek(x,50);
    m_est = Estymator_ortogonalny(X_N, Y_N, x, L);
    %m_est = Estymator_Fourier(X_N, Y_N, x, L, 7);
    suma = 0;
    for i=1:1:length(x)
       suma = suma + (real(m_est(i)) - m(i))^2; 
    end
    vErr = 1/length(x) * suma;
end

function m_est = Estymator_ortogonalny(X_N, Y_N, x, L)
    m_est = [];
    for i=1:1:length(x)
        pom = Estymator_f(X_N,L,x(i));
        if(pom == 0)
            m_est(end+1) = 0;
        else
            m_est(end+1) = Estymator_g(X_N,Y_N,L,x(i))/pom; 
        end
    end
end

function g_est = Estymator_g(X_N, Y_N, L, x)
    suma = 0;
    for k=0:1:L
        suma = suma + Estymator_alfa(X_N,Y_N,k)*Estymator_fi(k,x);
    end
    g_est = suma;
end

function alfa_est = Estymator_alfa(X_N,Y_N,k)
    suma = 0;
    for n=1:1:length(X_N)
        suma = suma + Y_N(n)* Estymator_fi(k,X_N(n));
    end
    alfa_est = 1/length(X_N) * suma;
end

function f_est = Estymator_f(X_N, L, x)
    suma = 0;
    for k=0:1:L
        suma = suma + Estymator_beta(X_N,k)*Estymator_fi(k,x);
    end
    f_est = suma;
end

function beta_est = Estymator_beta(X_N,k)
    suma = 0;
    for n=1:1:length(X_N)
       suma = suma + Estymator_fi(k,X_N(n)); 
    end
    beta_est = 1/length(X_N) * suma;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Bazy                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Baza cosinusowa
function fi = Estymator_fi(k,x)
    if(k == 0)
       fi = sqrt(1/(2*pi));
    else
       fi = sqrt(1/pi)*cos(k*x);
    end
end

%%%% Baza cosinusowa z elementami sinusa
% function fi = Estymator_fi(k,x)
%     if(k == 0)
%         fi = sqrt(1/(2*pi));
%     else
%         if(mod(k,2) == 0)
%             fi = sqrt(1/pi)*sin(k*x);
%         else
%             fi = sqrt(1/pi)*cos(k*x);
%         end
%     end
% end

%%%%% Baza Hermite'a
% function fi = Estymator_fi(k,x)
%     if(k == 0)
%         H_n = 1;
%     elseif(k == 1)
%         H_n = 2*x;
%     elseif(k == 2)
%         H_n = 4*x^2-2;
%     elseif(k == 3)
%         H_n = 8*x^3-12*x;
%     elseif(k == 4)
%         H_n = 16*x^4-48*x^2+12;
%     elseif(k == 5)
%         H_n = 32*x^5-160*x^3+120*x;
%     elseif(k == 6)
%         H_n = 64*x^6-480*x^4+720*x^2-120;
%     elseif(k == 7)
%         H_n = 128*x^7-1344*x^5+3360*x^3-1680*x;
%     end
%     fi = 1/sqrt(2^k*factorial(k)*sqrt(pi)) * exp(-x^2/2)* H_n;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Systemy                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y_N = Nieliniowy_system_statyczny(X_N, Z_N, a)
    f = @(x) (abs(x)>=2) .* (0) + (abs(x)>=0 && abs(x)<1) .* (a*x^2) + (abs(x)>=1 && abs(x)<2) .* (1);
    Y_N = zeros(1,length(X_N));
    for i=1:1:length(X_N)
        Y_N(i) = f(X_N(i)) + Z_N(i);
    end
end

function Y_N = System_statyczny_exp(X_N, Z_N)
    f = @(x) -x*exp(-x^2);
    Y_N = zeros(1,length(X_N));
    for i=1:1:length(X_N)
       Y_N(i) = f(X_N(i)) + Z_N(i); 
    end
end

function m = Zrob_funkcje_kotek(x, a)
    f = @(x) (abs(x)>=2) .* (0) + (abs(x)>=0 && abs(x)<1) .* (a*x^2) + (abs(x)>=1 && abs(x)<2) .* (1);
    m = [];
    for i=1:1:length(x)
        m(end+1) = f(x(i));
    end
end

function m = Zrob_funkcje_exp(x)
     f = @(x) -x*exp(-x^2);
     m = [];
    for i=1:1:length(x)
        m(end+1) = f(x(i));
    end
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

function setGlobal(x, y)
    global X_N;
    global Y_N;
    
    X_N = x;
    Y_N = y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Fourier                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function est = Estymator_Fourier(X_N, Y_N, x, L, T)
    est = [];
    for i=1:1:length(x)
        pom = Mianownik(X_N, x(i), L, T);
        if(pom == 0)
            est(end+1) = 0;
        else
            est(end+1) = Licznik(X_N,Y_N,x(i),L, T)/pom; 
        end
    end
end

function gora = Licznik(X_N, Y_N, x, L, T)
    suma = 0;
    omega = 2*pi/T;
    for k=-L:1:L
        suma = suma + Est_a(X_N, Y_N, k, T) * exp(1i*k*x*omega);
    end
    gora = suma;
end

function a = Est_a(X_N, Y_N, k, T)
    suma = 0;
    omega = 2*pi/T;
    for i=1:1:length(X_N)
        suma = suma + Y_N(i) * exp(-1i*k*X_N(i)*omega);
    end
    a = 1/(2*pi*length(X_N))*suma;
end

function dol = Mianownik(X_N, x, L, T)
    suma = 0;
    omega = 2*pi/T;
    for k=-L:1:L
        suma = suma + Est_b(X_N,k,T) * exp(1i*k*x*omega);
    end 
    dol = suma;
end

function b = Est_b(X_N,k, T)
    suma = 0;
    omega = 2*pi/T;
    for i=1:1:length(X_N)
        suma = suma + exp(-1i*k*X_N(i)*omega);
    end
    b = 1/(2*pi*length(X_N))*suma;
end

function Cvar = Cauchy(N)
    %cvar = @(N) tan((rand(1,N) - 0.5)*pi);
    cvar = @(N) randn(1,N)./randn(1,N);
    Cvar = cvar(N);
end
