%%
clear;
close all;
clc;

%% Stałe
s = 10;
N = 1000;
sigZ = 1;
b = ones(1,s)';
Un = normrnd(0,1,[1,N])';
Z = normrnd(0,sigZ,[1,N])';

%% Generowanie danych system SISO nieskorelowany
Zn = Z;
[Yn, PhiN] = System_SISO(Un, Zn, b);

%% Generowanie danych system SISO skorelowany
alfa = 0.5;
Zn = Zaklocenie_skor(Z, alfa);
[Yn, PhiN] = System_SISO(Un, Zn, b);

%% Estymator najmniejszych kwadratów
b_est = inv(PhiN' * PhiN)*PhiN'*Yn;
MSE = norm(b_est - b)^2;

figure(1);
plot(b, 'o');
hold on;
plot(b_est, '*');
title("Porownanie wektora parametrow systemu SISO z estymowanym metoda MNK - Err = " + MSE, 'interpreter','latex');
ylabel('b','interpreter','latex');

%% Blad empiryczny w funkcji N
L = 20;
Err = [];
for n=100:10:1000
    Err(end+1) = Blad_empiryczny(L, b, n, sigZ);
end
n = 100:10:1000;
figure(1);
plot(n, Err);
title('Blad empiryczny estymatora najmniejszych kwadratow w funkcji N','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');

%% Blad empiryczny z zakloceniem skorelowanym
L = 20;
alfa = 0.5;
Err_skor = [];
for n=100:10:1000
    Err_skor(end+1) = Blad_empiryczny_skor(L, b, n, sigZ, alfa);
end
n = 100:10:1000;
figure(1);
plot(n, Err_skor);
title('Blad empiryczny estymatora najmniejszych kwadratow w funkcji N','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');

%% Porownanie bledów 
figure(1);
plot(n, Err);
hold on;
plot(n, Err_skor);
title('Porownanie bledow identyfikacji w obecnosci skorelowanego i nieskorelowanego zaklocenia','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');
legend('nieskor','skor');

%% Macierz kowariancji zakłócenia
pom = zeros(1,N-2);
r = [(1+alfa^2)*sigZ, alfa*sigZ];
pom = [r,pom];
R = toeplitz(pom);
wyzn = det(R); % sprawdzenie czy osobliwa

imagesc(R);
colorbar;
title("Macierz kowariancji zaklocenia R, det(R) = " + wyzn,'interpreter','latex');

P = chol(R);

% Z_mat = [];
% for i=1:1:3000
%    Z_ = normrnd(0,sigZ,[1,N]);
%    Z_cor = inv(P) * Zaklocenie_skor(Z_, alfa);
%    Z_mat = [Z_mat, Z_cor];
% end
% C = cov(Z_mat');
% wyzn = det(C);
% imagesc(C);
% title("Macierz kowariancji zaklocenia R, det(R) = " + wyzn,'interpreter','latex');
% colorbar;

%% Wybielanie szumu?
Zn_ = inv(P) * Zn;

[Yn, PhiN] = System_SISO(Un, Zn_, b);
Yn_ =  Yn;
PhiN_ =  PhiN;

b_est = inv(PhiN_' * PhiN_)*PhiN_'*Yn_;
MSE = norm(b_est - b)^2;

figure(1);
plot(b, 'o');
hold on;
plot(b_est, '*');
title("Porownanie wektora parametrow systemu SISO z estymowanym metoda MNK - Err = " + MSE, 'interpreter','latex');
ylabel('b','interpreter','latex');

%% Błąd wybielonego estymatora

L = 20;
alfa = 0.5;
Err_wyb = [];
for n=100:10:1000
    Err_wyb(end+1) = Blad_empiryczny_wybielony(L, b, n, sigZ, alfa);
end
n = 100:10:1000;
figure(1);
plot(n, Err_wyb);
title('Blad empiryczny estymatora najmniejszych kwadratow w funkcji N','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');

%%
figure(1);

plot(n, Err_skor);
hold on;
plot(n, Err_wyb);

title('Porownanie bledow identyfikacji w obecnosci skorelowanego i wybielonego zaklocenia','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');
legend('skor','wyb');

%% Uogólniona metoda najmniejszych kwadratów
b_gls = inv(PhiN'*inv(R)*PhiN)*PhiN'*inv(R)*Yn;
MSE = norm(b_gls - b)^2;
figure(1);
plot(b, 'o');
hold on;
plot(b_gls, '*');
title("Porownanie wektora parametrow systemu SISO z estymowanym uogolniona metoda MNK - Err = " + MSE, 'interpreter','latex');
ylabel('b','interpreter','latex');

%% Blad empiryczny GLS
L = 20;
alfa = 0.5;
Err_gls = [];
for n=100:10:1000
    Err_gls(end+1) = Blad_empiryczny_gls(L, b, n, sigZ, alfa);
end
n = 100:10:1000;
figure(1);
plot(n, Err_gls);
title('Blad empiryczny estymatora najmniejszych kwadratow w funkcji N','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');

%% Porównanie wszystkich błędów
figure(1);
plot(n, Err);
hold on;
plot(n, Err_skor);
hold on;
plot(n, Err_gls);
legend('nieskor','skor', 'skor gls');
title('Porownanie bledow identyfikacji w obecnosci skorelowanego zaklocenia','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');

%% Zakłócenie z dodatkową korelacją
alfa = 0.5;
Zn = Zaklocenie_super(Z, alfa);
[Yn, PhiN] = System_SISO(Un, Zn, b);

%% Macierz kowariancji zakłócenia
Z_mat = [];
for i=1:1:3000
   Z_ = normrnd(0,sigZ,[1,N]);
   Z_cor = Zaklocenie_super(Z_, alfa);
   Z_mat = [Z_mat, Z_cor];
end
C = cov(Z_mat');
wyzn = det(C);
imagesc(C);
title("Macierz kowariancji zaklocenia R, det(R) = " + wyzn,'interpreter','latex');
colorbar;

%%
b_gls = inv(PhiN'*pinv(C)*PhiN)*PhiN'*pinv(C)*Yn;
MSE = norm(b_gls - b)^2;
figure(1);
plot(b, 'o');
hold on;
plot(b_gls, '*');
title("Porownanie wektora parametrow systemu SISO z estymowanym uogolniona metoda MNK - Err = " + MSE, 'interpreter','latex');
ylabel('b','interpreter','latex');

%% Blad empiryczny super zakłócenie
L = 20;
alfa = 0.5;
Err_sup = [];
for n=100:10:500
    Err_sup(end+1) = Blad_empiryczny_super(L, b, n, sigZ, alfa);
end

Err_sup_norm = [];
for n=100:10:500
    Err_sup_norm(end+1) = Blad_empiryczny_super_norm(L, b, n, sigZ, alfa);
end

n = 100:10:500;
figure(1);
plot(n, Err_sup_norm);
hold on;
plot(n, Err_sup);
legend('norm','gls');
title('Porownanie bledow identyfikacji w obecnosci skorelowanego zaklocenia','interpreter','latex');
xlabel('N','interpreter','latex');
ylabel('Err','interpreter','latex');

%% Funkcje małe i duże
function Z = Zaklocenie_super(Zn, alfa)
    Z = zeros(1, length(Zn))';
    for i=1:1:length(Zn)
       if(i == 1)
           Z(i) = 0;
       else
           Z(i) = alfa*Z(i-1) + Zn(i-1);
       end
    end
end

function Z = Zaklocenie_skor(Zn, alfa)
    Z = zeros(1, length(Zn))';
    for i=1:1:length(Zn)
       if(i == 1)
          Z(i) = Zn(i); 
       else
           Z(i) = Zn(i) + alfa*Zn(i-1);
       end
    end
end

function Err = Blad_empiryczny_wybielony(L, b, N, sig, alfa)
    sigZ = sig;
    suma = zeros(1,L);
    p = zeros(1,N-2);
    r = [(1+alfa^2)*sigZ, alfa*sigZ];
    p = [r,p];
    R = toeplitz(p);
    P = chol(R);
    pom = inv(P);
    for l=1:1:L
        Un = normrnd(0,sig,[1,N])';
        Z = normrnd(0,sig,[1,N])';
        Zn = pom * Zaklocenie_skor(Z, alfa);
        [Yn, PhiN] = System_SISO(Un, Zn, b);
        PhiN = pom*PhiN; %%%%%%%%%%%%%%%%%%%%%%%%%%%
        Yn = pom*Yn;
        b_est = inv(PhiN' * PhiN)*PhiN'*Yn;
        suma(l) = norm(b_est - b)^2;
    end
    Err = 1/L * sum(suma);
end

function Err = Blad_empiryczny_super(L, b, N, sig, alfa)
    suma = zeros(1,L);
    for l=1:1:L
        Un = normrnd(0,sig,[1,N])';
        Z = normrnd(0,sig,[1,N])';
        Zn = Zaklocenie_super(Z, alfa);
        [Yn, PhiN] = System_SISO(Un, Zn, b);
        
        Z_mat = [];
        for i=1:1:3000
           Z_ = normrnd(0,1,[1,N]);
           Z_cor = Zaklocenie_super(Z_, alfa);
           Z_mat = [Z_mat, Z_cor];
        end
        R = cov(Z_mat');
        
        b_est = inv(PhiN'*pinv(R)*PhiN)*PhiN'*pinv(R)*Yn;
        suma(l) = norm(b_est - b)^2;
    end
    Err = 1/L * sum(suma);
end

function Err = Blad_empiryczny_gls(L, b, N, sig, alfa)
    suma = zeros(1,L);
    for l=1:1:L
        Un = normrnd(0,sig,[1,N])';
        Z = normrnd(0,sig,[1,N])';
        Zn = Zaklocenie_skor(Z, alfa);
        [Yn, PhiN] = System_SISO(Un, Zn, b);
        
        pom = zeros(1,N-2);
        r = [(1+alfa^2)*sig, alfa*sig];
        pom = [r,pom];
        R = toeplitz(pom);
        
        b_est = inv(PhiN'*inv(R)*PhiN)*PhiN'*inv(R)*Yn;
        suma(l) = norm(b_est - b)^2;
    end
    Err = 1/L * sum(suma);
end

function Err = Blad_empiryczny_super_norm(L, b ,N, sig, alfa)
    suma = zeros(1,L);
    for l=1:1:L
        Un = normrnd(0,sig,[1,N])';
        Z = normrnd(0,sig,[1,N])';
        Zn = Zaklocenie_super(Z, alfa);
        [Yn, PhiN] = System_SISO(Un, Zn, b);
        b_est = inv(PhiN' * PhiN)*PhiN'*Yn;
        suma(l) = norm(b_est - b)^2;
    end
    Err = 1/L * sum(suma);
end

function Err = Blad_empiryczny_skor(L, b, N, sig, alfa)
    suma = zeros(1,L);
    for l=1:1:L
        Un = normrnd(0,sig,[1,N])';
        Z = normrnd(0,sig,[1,N])';
        Zn = Zaklocenie_skor(Z, alfa);
        [Yn, PhiN] = System_SISO(Un, Zn, b);
        b_est = inv(PhiN' * PhiN)*PhiN'*Yn;
        suma(l) = norm(b_est - b)^2;
    end
    Err = 1/L * sum(suma);
end

function Err = Blad_empiryczny(L, b, N, sig)
    suma = zeros(1,L);
    for l=1:1:L
        Un = normrnd(0,sig,[1,N])';
        Zn = normrnd(0,sig,[1,N])';
        [Yn, PhiN] = System_SISO(Un, Zn, b);
        b_est = inv(PhiN' * PhiN)*PhiN'*Yn;
        suma(l) = norm(b_est - b)^2;
    end
    Err = 1/L * sum(suma);
end

function [Yn, PhiN] = System_SISO(Un, Zn, b)
    phi = zeros(1,length(b));
    PhiN = [];
    for i=1:1:length(Un)
        phi = [Un(i), phi];
        phi(end) = [];
        PhiN = [PhiN; phi];
    end
    Yn = PhiN * b + Zn;
end
