%%
clear;
close all;
clc;

%% Definicja systemu i generacja sygnałów
N = 100; % - ilość danych pomiarowych
D = 120; % - ilość wejść systemu
a = ones(1,D)'; % - wektor a opisujący system
%a = rand(1,D)'; % - wektor a opisujący system


% Generacja sygnałów wejściowych
sigX = 4;
Sig = eye(D) * sigX;
mu = rand(1,D);
%mu = 2*ones(1,D);
Xn = mvnrnd(mu,Sig,N);

% Generacja szumu i wyjścia systemu
sigZ = 1;
Zn = normrnd(0,sigZ,[1,N])';
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

%% Graficzna reprezentacja macierzy kowariancji
Cov = sigZ*inv(Xn' * Xn);
imagesc(Cov);
colorbar;
title("Graficzna reprezentacja macierzy kowariancji N=100, D=100",'interpreter','latex');

%% Błąd empiryczny
L = 20;

suma = 0;
Err = [];

for n=100:10:1000
    suma = 0;
    for l=1:1:L
        D = round(n/3);
        a = ones(1,D)'; % - wektor a opisujący system

        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = 2*ones(1,D);
        Xn = mvnrnd(mu,Sig,n);

        % Generacja szumu i wyjścia systemu
        sigZ = 1;
        Zn = normrnd(0,sigZ,[1,n])';
        Yn = System_MISO(Xn, Zn, a);

        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;

        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err(end+1) = 1/L * suma;
end

n_ = 100:10:1000;
figure(1);
plot(n_, Err);
title('Blad empiryczny w funkcji N dla roznych wariancji zaklocenia','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');
legend('sig=1','sig=5','sig=10','interpreter','latex')

%%
N = 1000; % - ilość danych pomiarowych
D = N; % - ilość wejść systemu
a = ones(1,D)'; % - wektor a opisujący system
%a = rand(1,D)'; % - wektor a opisujący system

% Generacja sygnałów wejściowych
sigX = 4;
Sig = eye(D) * sigX;
mu = rand(1,D);
%mu = 2*ones(1,D);
Xn = mvnrnd(mu,Sig,N);

% Generacja szumu i wyjścia systemu
sigZ = 1;
Zn = normrnd(0,sigZ,[1,N])';
Yn = System_MISO(Xn, Zn, a);

a_est = inv(Xn' * Xn)*Xn'*Yn;
MSE = norm(a_est - a)^2;

L = 20;
Err = [];
for n = 100:50:1000
    Err(end+1) = Blad_empiryczny(L, n);
end
n_ = 100:50:1000;

subplot(2,1,1);
figure(1);
plot(a,'*');
hold on;
plot(a_est,'o');
title("Porownanie wektora parametrow systemu MISO z estymowanym metoda MNK - Err = " + MSE, 'interpreter','latex');
ylabel('a','interpreter','latex');
subplot(2,1,2)
plot(n_,Err,'*');
title('Blad empiryczny estymatora najmniejszych kwadratow w funkcji N - D = N','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');

%% Estymacja macierzy kowariancji estymatora

% Estymacja wariancji zakłócenia
sig_est = 1/(N - D)*norm(Yn - Xn*a_est)^2;

% Estymacja macierzy kowariancji
Cov_est = sig_est * inv(Xn'*Xn);

% figure;
% imagesc(Cov_est);
% colorbar;
% title("Macierz kowariancji zbudowana na podstawie estymacji wariancji zaklocenia sig: " + sig_est,'interpreter','latex');
diff = abs(Cov - Cov_est);
imagesc(diff);
colorbar;
title('Bledy estymacji macierzy kowariancji','interpreter','latex');

%% Estymacja macierzy kowariancji - dane wejściowe co się dzieje gdy są zależne od siebie wielkości
sum = 0;
for i=1:1:length(Yn)
    sum = sum + ( (Xn(i,:) - mean(Xn(i,:)))'*(Xn(i,:) - mean(Xn(i,:))) );
end
Cov_est = 1/length(Yn) * sum;

%% Funkcje
function Err = Blad_empiryczny(L, N)
    suma = 0;
    for l=1:1:L
        
        D = round(N);
        %D = 10; % - ilość wejść systemu
        a = ones(1,D)'; % - wektor a opisujący system
        
        % Generacja sygnałów wejściowych
        sigX = 4;
        Sig = eye(D) * sigX;
        mu = 2*ones(1,D);
        Xn = mvnrnd(mu,Sig,N);

        % Generacja szumu i wyjścia systemu
        sigZ = 1;
        Zn = normrnd(0,sigZ,[1,N])';
        Yn = System_MISO(Xn, Zn, a);
        
        % Estymacja
        a_est = inv(Xn' * Xn)*Xn'*Yn;
        
        % Wyznaczenie błędu realizacji
        suma = suma + norm(a_est - a)^2;
    end
    Err = 1/L * suma;
end

function Yn = System_MISO(Xn, Zn, a)
    Yn = zeros(1,length(Zn));
    for i=1:1:length(Zn)
        Yn(i) = Xn(i,:) * a;
    end
    Yn = Yn' + Zn;
end
