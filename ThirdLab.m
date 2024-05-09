close all;
clear all;
coef = 1; #регулировка четкости картинки
lambda = 0.532 * coef; % длина волны
k = (2 * pi)/lambda; % волновое число
window_size = 4; % размер получаемых изображений
R = 2.5 * lambda;
z1 = 0.5*lambda; % расстояния от экрана
z2 = lambda;
z3 = 5*lambda;
z4 = 10*lambda;
z5 = 50*lambda;

% матрицы для результирующих полей
U1 = zeros(2*window_size*coef, 2*window_size*coef);
U2 = zeros(2*window_size*coef, 2*window_size*coef);
U3 = zeros(2*window_size*coef, 2*window_size*coef);
U4 = zeros(2*window_size*coef, 2*window_size*coef);


for i = 1:(2*window_size*coef + 1)
  for j = 1:(2*window_size*coef + 1)

    X = i-(window_size*coef + 1);
    Y = j-(window_size*coef + 1);

    % подынтегральные функции
    int1 = @(x,y) exp(1i*k*z1).*(exp(1i.*k.*(sqrt((X-x).^2+(Y-y).^2+z1.^2)))./((X-x).^2+(Y-y).^2+z1.^2)).*(1+(1i./(k.*sqrt((X-x).^2+(Y-y).^2+z1.^2))));
    int2 = @(x,y) exp(1i*k*z2).*(exp(1i.*k.*(sqrt((X-x).^2+(Y-y).^2+z2.^2)))./((X-x).^2+(Y-y).^2+z2.^2)).*(1+(1i./(k.*sqrt((X-x).^2+(Y-y).^2+z2.^2))));
    int3 = @(x,y) exp(1i*k*z3).*(exp(1i.*k.*(sqrt((X-x).^2+(Y-y).^2+z3.^2)))./((X-x).^2+(Y-y).^2+z3.^2)).*(1+(1i./(k.*sqrt((X-x).^2+(Y-y).^2+z3.^2))));
    int4 = @(x,y) exp(1i*k*z4).*(exp(1i.*k.*(sqrt((X-x).^2+(Y-y).^2+z4.^2)))./((X-x).^2+(Y-y).^2+z4.^2)).*(1+(1i./(k.*sqrt((X-x).^2+(Y-y).^2+z4.^2))));
    int5 = @(x,y) exp(1i*k*z5).*(exp(1i.*k.*(sqrt((X-x).^2+(Y-y).^2+z5.^2)))./((X-x).^2+(Y-y).^2+z5.^2)).*(1+(1i./(k.*sqrt((X-x).^2+(Y-y).^2+z5.^2))));
    % вычисление интегралов Рэлея-Зоммерфельда первого рода
    U1(i,j) = integral2(int1,-R,R, @(j)-sqrt(R^2 - j.^2), @(j)sqrt(R^2 - j.^2));%круглая апертура
    U2(i,j) = integral2(int2,-R,R, @(j)-sqrt(R^2 - j.^2), @(j)sqrt(R^2 - j.^2));
    U3(i,j) = integral2(int3,-R,R, @(j)-sqrt(R^2 - j.^2), @(j)sqrt(R^2 - j.^2));
    U4(i,j) = integral2(int4,-R,R, @(j)-sqrt(R^2 - j.^2), @(j)sqrt(R^2 - j.^2));
    U5(i,j) = integral2(int5,-R,R, @(j)-sqrt(R^2 - j.^2), @(j)sqrt(R^2 - j.^2));
  end
end

U1 = (-(1i*k*z1)/(2*pi)).*U1;
U2 = (-(1i*k*z2)/(2*pi)).*U2;
U3 = (-(1i*k*z3)/(2*pi)).*U3;
U4 = (-(1i*k*z4)/(2*pi)).*U4;
U5 = (-(1i*k*z5)/(2*pi)).*U5;
#вычисляем интенсивность
center_row = round(size(U1, 1) / 2);
intensity = U1(center_row, :) .* conj(U1(center_row, :));
x_axis = linspace(-window_size, window_size, length(intensity));


figure(1);
plot(x_axis, real(intensity));
xlabel('x');
ylabel('Интенсивность');
title('z = 0.5*lambda');


center_row = round(size(U2, 1) / 2);
intensity = U2(center_row, :) .* conj(U2(center_row, :));
x_axis = linspace(-window_size, window_size, length(intensity));

figure(2);
plot(x_axis, real(intensity));
xlabel('x');
ylabel('Интенсивность');
title('z = lambda');


center_row = round(size(U3, 1) / 2);
intensity = U3(center_row, :) .* conj(U3(center_row, :));
x_axis = linspace(-window_size, window_size, length(intensity));

figure(3);
plot(x_axis, real(intensity));
xlabel('x');
ylabel('Интенсивность');
title('z = 5*lambda');


center_row = round(size(U4, 1) / 2);
intensity = U4(center_row, :) .* conj(U4(center_row, :));
x_axis = linspace(-window_size, window_size, length(intensity));

figure(4);
plot(x_axis, real(intensity));
xlabel('x');
ylabel('Интенсивность');
title('z = 10*lambda');


center_row = round(size(U5, 1) / 2);
intensity = U5(center_row, :) .* conj(U5(center_row, :));
x_axis = linspace(-window_size, window_size, length(intensity));

figure(5);
plot(x_axis, real(intensity));
xlabel('x');
ylabel('Интенсивность');
title('z = 50*lambda');

%figure(6);
%plot(r, airy_10);
%xlabel('x');
%ylabel('Интенсивность');
%title('z = 10*lambda (airy)');


%figure(7);
%plot(r, airy_50);
%xlabel('x');
%ylabel('Интенсивность');
%title('z = 50*lambda (airy)');



r = linspace(-window_size, window_size, 100)#распределение
%Дифракция Фраунгофера
airy_10 = (besselj(1, k * 2 * R * r / (2 * z4)) ./ (k * 2 * R * r / (2 * z4))).^2;
airy_50 = (besselj(1, k * 2 * R * r / (2 * z5)) ./ (k * 2 * R * r / (2 * z5))).^2;

% Обработка случая, когда знаменатель равен нулю (в центре)
airy_10(r == 0) = 1;
airy_50(r == 0) = 1;

center_row = round(size(U4, 1) / 2);
intensity = U4(center_row, :) .* conj(U4(center_row, :));
r = linspace(0, 3, length(airy_10));
original_x = linspace(0, 3, length(intensity));
interp_intensity = interp1(original_x, intensity, r, 'linear', 'extrap');
normalized_airy = airy_10 / max(airy_10);
normalized_intensity = interp_intensity / max(interp_intensity);


figure(6);
plot(r, normalized_intensity, 'b', 'DisplayName', 'RS-I');
hold on;
plot(r, normalized_airy, 'r', 'DisplayName', 'Airy');
xlabel('r');
ylabel('Интенсивность');
title('Сравнение интенсивности (RS-1, Airy) при z = 10lambda');
legend show;
% Настройка границ осей
%axis([0, 3, 0, 1]); % [xmin xmax ymin ymax]


center_row = round(size(U5, 1) / 2);
intensity = U5(center_row, :) .* conj(U5(center_row, :));
r = linspace(0, 15, length(airy_50));
original_x = linspace(0, 15, length(intensity));
interp_intensity = interp1(original_x, intensity, r, 'linear', 'extrap');
normalized_airy = airy_50 / max(airy_50);
normalized_intensity = interp_intensity / max(interp_intensity);


figure(7);
plot(r, normalized_intensity, 'b', 'DisplayName', 'RS-I');
hold on;
plot(r, normalized_airy, 'r', 'DisplayName', 'Airy');
xlabel('r');
ylabel('Интенсивность');
title('Сравнение интенсивности (RS-1, Airy) при z = 50lambda');
legend show;
grid on;
axis([0, 15, 0, 1]);


x=-window_size:window_size;
y=-window_size:window_size;
figure(8), imagesc(x,y,U1.*conj(U1)); colorbar; xlabel('\it{x, \mum}'); ylabel('\it{y, \mum}'); set(gca,'YDir','normal'); title('\it{z = \0.5lambda}');
