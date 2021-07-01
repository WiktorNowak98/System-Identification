%% 
clear;
close all;
clc;

%% Stałe
N = 500;
a = 0.8; % a musi być z przedziału -1,1 aby system był stabilny ujemne wartości wprowadzają oscylacje
b = 5;
theta = [b, a]';

%% Zabawa odpowiedzią impulsową systemu dynamicznego
Imp = [1, zeros(1, N-1)];
Res = 0;
Yn = [];
for i=1:1:N
    fi = [Imp(i), Res]';
    Res = fi' * theta;
    Yn(end+1) = Res;
end

figure(1);
plot(0:N-1, Yn);
xlabel('Xn');
ylabel('Yn');
title("Odpowiedź impulsowa systemu - parametry a = " + a + " b = " + b);

%% Generacja pomiarów wejścia wyjścia systemu dynamicznego
Un = rand(1, N);
en = 2*rand(1,N)-1;

[Yn, Phi] = System(Un, en, theta, a);

%% Estymator najmniejszych kwadratów
theta_est = inv(Phi' * Phi) * Phi' * Yn;
MSE = norm(theta_est - theta)^2;
figure(1);
plot(theta,'o');
hold on;
plot(theta_est,'x');
title("Przykładowa realizacja estymatora najmniejszych kwadratów MSE = " + MSE);
    
%% Błąd empiryczny estymatora najmniejszych kwadratów
L = 50;
Err_MNK = [];
for n = 10:1:500
    Err_MNK(end+1) = EmpErrMNK(L, n, theta, a);
end
n = 10:1:500;
plot(n, Err_MNK);
title('Błąd empiryczny estymatora najmnhiejszych kwadratów w funkcji N')

%% Oszacowanie niezaszumionego wyjścia V
Vn_est = [];
pom = 0;
for i=1:1:N
    Vn_est(end+1) = theta_est(1) * Un(i) + theta_est(2) * pom;
    pom = Vn_est(end);
end
Vn_est = [0,Vn_est];
Vn_est(end) = [];

Psi = [Un',Vn_est'];

%% Metoda zmiennych instrumentalnych
theta_IV = inv(Psi' * Phi)* Psi' * Yn;
MSE = norm(theta_IV - theta)^2;
figure(1);
plot(theta,'o');
hold on;
plot(theta_IV,'x');
title("Przykładowa realizacja estymatora zmiennych instrumentalnych MSE = " + MSE);
%% Błąd empiryczny zmiennych instrumentalnych
L = 50;
Err_IV = [];
for n = 10:1:500
    Err_IV(end+1) = EmpErrIV(L, n, theta, a);
end
n = 10:1:500;
plot(n, Err_IV);

%% Porównanie błędów
figure(1);
plot(n, Err_MNK);
hold on;
plot(n, Err_IV);
legend('MNK','IV');
title('Porównanie błędu empirycznego metody najmniejszych kwadratów i metody zmiennych instrumentalnych');
%% Funkcje małe i duże

function [Yn, Phi] = System(Un, en, theta, a)
    Zn = [];
    Res = 0;
    Phi = [];
    e_s = 0;

    for i=1:1:length(Un)
        fi = [Un(i), Res]';
        Phi = [Phi;fi'];
        Res = fi' * theta + en(i);
        Zn(end+1) = en(i) - a*e_s;
        e_s = en(i);
    end
    Yn = Phi * theta + Zn';
end

function Err = EmpErrMNK(L, N, theta, a)
    suma = 0;
    for l=1:1:L
       Un = rand(1, N);
       en = 2*rand(1,N)-1;
       [Yn, Phi] = System(Un, en, theta, a);
       theta_est = inv(Phi' * Phi) * Phi' * Yn;
       suma = suma + norm(theta_est - theta)^2;
    end
    Err = 1/L * suma;
end

function Err = EmpErrIV(L, N, theta, a)
    suma = 0;
    for l=1:1:L
       Un = rand(1, N);
       en = 2*rand(1,N)-1;
       [Yn, Phi] = System(Un, en, theta, a);
       theta_est = inv(Phi' * Phi) * Phi' * Yn;
       
       Vn_est = [];
       pom = 0;
       for i=1:1:N
           Vn_est(end+1) = theta_est(1) * Un(i) + theta_est(2) * pom;
           pom = Vn_est(end);
       end
       Vn_est = [0,Vn_est];
       Vn_est(end) = [];
       Psi = [Un',Vn_est'];
       
       theta_IV = inv(Psi' * Phi)* Psi' * Yn;
       
       suma = suma + norm(theta_IV - theta)^2;
    end
    Err = 1/L * suma;
end





