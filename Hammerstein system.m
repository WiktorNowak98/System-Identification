%%
clear;
close all;
clc;

%% Dane
N = 1000;
L = 25;
b = .8*ones(1,L);
for l=1:1:L
    b(l) = b(l)^l;
end
b = b';
Un = 2*rand(1, N)-1;
Zn = normrnd(0,1,[1,N]);
Yn = Hammerstein(Un, Zn, b)';

x = -1:0.01:1;
y = Nielinowy_statyczny(x);

%% Estymator nieliniowości systemu Hammersteina
jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
h = 0.5;
m_est = Estymator_jadrowy(Un, Yn, h, jadro_box, x);
mse = mean(m_est - y)^2;

subplot(2,1,1);
plot(Un, Yn, '.');
hold on;
plot(x, m_est);
title('Chmura danych pomiarowych wraz z wynikiem estymacji');
legend('dane','est');
subplot(2,1,2);
plot(x, y);
hold on;
plot(x, m_est);
title("Estymacja nieliniowości systemu Hammersteina MSE = " + mse);
legend('sys', 'est');

%% Błąd empiryczny nieliniowej
Err_ = [];
for n = 10:10:1000
    Err_(end+1) = Blad_niel(n, 20);
end
n = 10:10:1000;
figure(1);
plot(n, Err_);
title("Błąd empiryczny estymacji nieliniowości systemu Hammersteina w funkcji N");

%% Estymator podsystemu liniowego
b_est = Estymator_dynamiki(Un, Yn, b)';
mse_l = norm(b_est - b)^2;

figure(1);
plot(b,'o');
hold on;
plot(b_est, '*');
title("Estymacja odpowiedzi impulsowej dynamicznego podsystemu liniowego");
legend('org','est');

%% Błąd empiryczny liniowej
Err = [];
for n = 10:50:2500
    Err(end+1) = Blad_empiryczny(n, 20, b);
end
n = 10:50:2500;
figure(1);
plot(n, Err);
title("Błąd empiryczny estymacji dynamicznego podsystemu liniowego");

%% Funkcje

function Err = Blad_niel(N, L)
    x = -1:0.01:1;
    y = Nielinowy_statyczny(x);
    b = ones(1, 5)';
    jadro_box = @(x) (abs(x) <= 0.5) .* (1) + (abs(x) > 0.5) .* (0);
    h = 0.5;
    suma = 0;
    for l=1:1:L
        Un = 2*rand(1, N)-1;
        Zn = normrnd(0,1,[1,N]);
        Yn = Hammerstein(Un, Zn, b)';
        m_est = Estymator_jadrowy(Un, Yn, h, jadro_box, x);
        suma = suma + mean(m_est - y)^2;
    end
    Err = 1/L * suma;
end

function Err = Blad_empiryczny(N, L, b)
    suma = 0;
    for l=1:1:L
        Un = 2*rand(1, N)-1;
        Zn = normrnd(0,1,[1,N]);
        Yn = Hammerstein(Un, Zn, b)';
        b_est = Estymator_dynamiki(Un, Yn, b)';
        suma = suma + norm(b_est - b)^2;
    end
    Err = 1/L * suma;
end

function b_est = Estymator_dynamiki(Un, Yn, b)
    b_est = [];
    for l=1:1:length(b)
       suma = 0;
       for k = l+1:1:length(Un)
           suma = suma + Un(k-l) * Yn(k);
       end
       b_est(end+1) = 1/(length(Un) - l) * suma;
    end
end

function m_est = Estymator_jadrowy(X_N, Y_N, h, jadro, x)
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

function Yn = Hammerstein(Un, Zn, b)
    Wn = Nielinowy_statyczny(Un);
    Yn = Liniowy_dynamiczny(Wn', Zn', b);
end

function Wn = Nielinowy_statyczny(Un)
    Wn = atan(12*Un);
end

function Yn = Liniowy_dynamiczny(Wn, Zn, b)
    phi = zeros(1,length(b));
    PhiN = [];
    for i=1:1:length(Wn)
        phi = [Wn(i), phi];
        phi(end) = [];
        PhiN = [PhiN; phi];
    end
    Yn = PhiN * b + Zn;
end