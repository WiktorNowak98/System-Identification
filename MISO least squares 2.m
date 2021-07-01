%%
clear;
close all;
clc;

%% Generowanie danych
N = 250;
D = 10;
a = ones(1,D)';

sigX = 4;
Sig = eye(D) * sigX;
mu = rand(1,D);
%mu = rand(1,D);
Xn = mvnrnd(mu,Sig,N);

sigZ = 1;
Z = normrnd(0,sigZ,[1,N]);

b = 0.5;
Zn = Zaklocenie(Z, b, N)';
Yn = System_MISO(Xn, Zn, a);

%% Estymator najmniejszych kwadratów
a_est = inv(Xn' * Xn)*Xn'*Yn;
MSE = norm(a_est - a)^2;

figure(1);
plot(a,'*');
hold on;
plot(a_est,'o');
title("Porownanie wektora parametrow systemu MISO z estymowanym metoda MNK - Err = " + MSE, 'interpreter','latex');
ylabel('a','interpreter','latex');

%% Macierz kowariancji estymatora
% Macierz kowariancji zakłócenia
pom = zeros(1,N-2);
r = [(1+b^2)*sigZ, b*sigZ];
pom = [r,pom];
R = toeplitz(pom);

figure(1);
imagesc(R);

Cov = inv(Xn'*Xn)*Xn'*R*Xn*inv(Xn'*Xn);

Z_mat = [];
for i=1:1:3000
   Z = normrnd(0,sigZ,[1,N]);
   Z_cor = Zaklocenie(Z, b, N)';
   Z_mat = [Z_mat, Z_cor];
end

C = cov(Z_mat');

% subplot(2,1,1);
% imagesc(R);
% title("Analitycznie wyznaczona macierz kowariancji skorelowanego zaklocenia", 'interpreter','latex');
% colorbar;
% subplot(2,1,2);
% imagesc(C);
% title("Symulacyjnie wyznaczona macierz kowariancji skorelowanego zaklocenia", 'interpreter','latex');
% colorbar;

imagesc(Cov);
title("Graficzna reprezentacja macierzy kowariancji estymatora",'interpreter','latex');
colorbar;

%% Błąd empiryczny w funkcji N
L = 20;
Err = [];
for n = 100:50:1000
    Err(end+1) = Blad_empiryczny(L, n);
end
n_ = 100:50:1000;
figure(1);
plot(n_,Err);
title('Blad empiryczny estymatora najmniejszych kwadratow w funkcji N','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');

%% Empiryczny znowu

L = 20;

suma = 0;
Err = [];
Err1 = [];
Err2 = [];

for n=100:10:1000
    suma = 0;
    for l=1:1:L
        D = 10; % - ilość wejść systemu
        a = ones(1,D)'; % - wektor a opisujący system
        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = rand(1,D);
        Xn = mvnrnd(mu,Sig,n);
        % Generacja szumu i wyjścia systemu
        sigZ = 1;
        b = 3;
        Z = normrnd(0,sigZ,[1,n])';
        Zn = Zaklocenie(Z, b, n)';
        Yn = System_MISO(Xn, Zn, a);
        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;
        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err(end+1) = 1/L * suma;
end

for n=100:10:1000
    suma = 0;
    for l=1:1:L
        D = 10; % - ilość wejść systemu
        a = ones(1,D)'; % - wektor a opisujący system
        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = rand(1,D);
        Xn = mvnrnd(mu,Sig,n);
        % Generacja szumu i wyjścia systemu
        sigZ = 5;
        b = 3;
        Z = normrnd(0,sigZ,[1,n])';
        Zn = Zaklocenie(Z, b, n)';
        Yn = System_MISO(Xn, Zn, a);
        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;
        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err1(end+1) = 1/L * suma;
end

for n=100:10:1000
    suma = 0;
    for l=1:1:L
        D = 10; % - ilość wejść systemu
        a = ones(1,D)'; % - wektor a opisujący system
        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = rand(1,D);
        Xn = mvnrnd(mu,Sig,n);
        % Generacja szumu i wyjścia systemu
        sigZ = 10;
        b = 3;
        Z = normrnd(0,sigZ,[1,n])';
        Zn = Zaklocenie(Z, b, n)';
        Yn = System_MISO(Xn, Zn, a);
        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;
        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err2(end+1) = 1/L * suma;
end

n_ = 100:10:1000;
figure(1);
plot(n_, Err);
hold on;
plot(n_, Err1);
hold on;
plot(n_, Err2);
title('Blad empiryczny w funkcji N dla roznych wariancji zaklocenia','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');
legend('sig=1','sig=5','sig=10','interpreter','latex')

%% Dodatkowe
ff = 0.98;
sigZ = 1;
N = 1000;
Z_mat = [];
K = 500;
lambda=[];
for k=0:1:K
    lambda(end+1) = ff^k;
end

for i=1:1:3000
   Z = normrnd(0,sigZ,[1,N]);
   Z_cor = Zaklocenie_ff(Z, ff, N, K)';
   Z_mat = [Z_mat, Z_cor];
end

% C = cov(Z_mat');
% %imagesc(C);
% bar3(C);
% title("Macierz kowariancji skorelowanego zaklocenia ff: " + ff + " K: " + K, 'interpreter','latex');
% colorbar;

%% Dodatkowe 2

clear;
close all;
clc;

N = 250;
D = 10;
a = ones(1,D)';

sigX = 4;
Sig = eye(D) * sigX;
mu = rand(1,D);
Xn = mvnrnd(mu,Sig,N);

ff = 0.5;
sigZ = 1;
Z_mat = [];
K = 15;
L = 20;
Z = normrnd(0,sigZ,[1,N])';
Zn = Zaklocenie_ff(Z, ff, N, K)';
Yn = System_MISO(Xn, Zn, a);
a_est = inv(Xn' * Xn)*Xn'*Yn;
MSE = norm(a_est - a)^2;

for i=1:1:3000
   Z = normrnd(0,sigZ,[1,N]);
   Z_cor = Zaklocenie_ff(Z, ff, N, K)';
   Z_mat = [Z_mat, Z_cor];
end
R = cov(Z_mat');
Cov = inv(Xn'*Xn)*Xn'*R*Xn*inv(Xn'*Xn);

Err = [];
for n=100:10:1000
    suma = 0;
    for l=1:1:L
        D = 10; % - ilość wejść systemu
        a = ones(1,D)'; % - wektor a opisujący system
        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = rand(1,D);
        Xn = mvnrnd(mu,Sig,n);
        % Generacja szumu i wyjścia systemu
        sigZ = 1;
        b = 0.5;
        Z = normrnd(0,sigZ,[1,n])';
        Zn = Zaklocenie_ff(Z, ff, n, K)';
        Yn = System_MISO(Xn, Zn, a);
        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;
        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err(end+1) = 1/L * suma;
end

Err1 = [];
for n=100:10:1000
    suma = 0;
    for l=1:1:L
        D = 10; % - ilość wejść systemu
        a = ones(1,D)'; % - wektor a opisujący system
        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = rand(1,D);
        Xn = mvnrnd(mu,Sig,n);
        % Generacja szumu i wyjścia systemu
        sigZ = 5;
        b = 0.5;
        Z = normrnd(0,sigZ,[1,n])';
        Zn = Zaklocenie_ff(Z, ff, n, K)';
        Yn = System_MISO(Xn, Zn, a);
        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;
        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err1(end+1) = 1/L * suma;
end

Err2 = [];
for n=100:10:1000
    suma = 0;
    for l=1:1:L
        D = 10; % - ilość wejść systemu
        a = ones(1,D)'; % - wektor a opisujący system
        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = rand(1,D);
        Xn = mvnrnd(mu,Sig,n);
        % Generacja szumu i wyjścia systemu
        sigZ = 10;
        b = 0.5;
        Z = normrnd(0,sigZ,[1,n])';
        Zn = Zaklocenie_ff(Z, ff, n, K)';
        Yn = System_MISO(Xn, Zn, a);
        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;
        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err2(end+1) = 1/L * suma;
end


n_ = 100:10:1000;

subplot(3,1,1);
plot(a,'*');
hold on;
plot(a_est,'o');
title("Porownanie wektora parametrow systemu MISO z estymowanym metoda MNK - Err = " + MSE, 'interpreter','latex');
ylabel('a','interpreter','latex');
subplot(3,1,2);
imagesc(Cov);
title("Graficzna reprezentacja macierzy kowariancji estymatora",'interpreter','latex');
colorbar;
subplot(3,1,3);
plot(n_, Err);
hold on;
plot(n_, Err1);
hold on;
plot(n_, Err2);
title('Blad empiryczny w funkcji N','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');
legend('sig=1','sig=5','sig=10','interpreter','latex')
%% Funkcje

function Err = Blad_empiryczny(L, N)
    suma = 0;
    for l=1:1:L
        D = 10; % - ilość wejść systemu
        a = ones(1,D)'; % - wektor a opisujący system
        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = rand(1,D);
        Xn = mvnrnd(mu,Sig,N);
        % Generacja szumu i wyjścia systemu
        sigZ = 1;
        b = 0.5;
        Z = normrnd(0,sigZ,[1,N])';
        Zn = Zaklocenie(Z, b, N)';
        Yn = System_MISO(Xn, Zn, a);
        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;
        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err = 1/L * suma;
end

function Zn = Zaklocenie_ff(Z, ff, N, K)
    Zn = zeros(1,N);
    lambda = [];
    for k=0:1:K
        lambda(end+1) = ff^k;
    end
    pom = zeros(1,length(lambda));
    for i=1:1:N
        pom = [Z(i),pom];
        pom(end) = [];
        Zn(i) = sum(lambda .* pom);
    end
end

function Zn = Zaklocenie(Z, b, N)
    Zn = zeros(1,N);
    for i=1:1:N
        if(i == 1)
           Zn(i) = Z(i);
        else
            Zn(i) = Z(i) + b*Z(i-1);
        end
    end
end

function Zn = Zaklocenie_skor(Z, b, N)
    Zn = zeros(1,N);
    pom = zeros(1,length(b));
    for i=1:1:N
        pom = [Z(i),pom];
        pom(end) = [];
        Zn(i) = sum(b .* pom);
    end
end

function Yn = System_MISO(Xn, Zn, a)
    Yn = zeros(1,length(Zn));
    for i=1:1:length(Zn)
        Yn(i) = Xn(i,:) * a;
    end
    Yn = Yn' + Zn;
end
