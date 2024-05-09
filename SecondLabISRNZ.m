close all;
clear all;
% Параметры задачи
lambda = 0.532; % Длина волны
n = 1; % Показатель преломления среды
k = 2*pi/lambda; % Волновое число


for NA_val = [0.65, 0.8, 0.95]

  alpha = asin(NA_val) %максимальный угол, определяемый числовой апертурой
  T = @(theta) cos(theta).^0.5; %для апланатического объектива
  r = 0;
  z = linspace(-3, 3, 100);
  psi = 0;

  %для плоской линейнополяризованной волны
  Ex_1 = @(r, z, theta) 1 * T(theta) * sin(theta) * (1 - cos(theta)) * exp(1i * k * z * cos(theta)) * besselj(2, k .* r .* sin(theta));
  Ex_2 = @(r, z, theta) 1 * T(theta) * sin(theta) * (1 + cos(theta)) * exp(1i * k * z * cos(theta)) * besselj(0, k .* r .* sin(theta));
  E_x = -1i * cos(2 * psi) * integral(@(theta) Ex_1(r, z, theta), 0, alpha, 'ArrayValued', true) - ...
  1i * integral(@(theta) Ex_2(r, z, theta), 0, alpha, 'ArrayValued', true) %вычисляет интеграл для каждого элемента массива, возвращает массив значений

  Ey_1 = @(r, z, theta) 1 * T(theta) * sin(theta) * (1 - cos(theta)) * exp(1i * k * z * cos(theta)) * besselj(2, k .* r .* sin(theta));
  E_y = -1i * sin(2*psi) * integral(@(theta) Ey_1(r, z, theta), 0, alpha, 'ArrayValued', true)

  Ez_1 = @(r, z, theta) 1 * T(theta) * sin(theta)^2 * exp(1i * k * z * cos(theta)) * besselj(1, k .* r .* sin(theta));
  E_z = -2 * cos(psi) * integral(@(theta) Ez_1(r, z, theta), 0, alpha, 'ArrayValued', true)

  I = abs(E_x).^2 + abs(E_y).^2 + abs(E_z).^2 %интенсивность считаем как сумму квадратов модулей

  % Визуализация результатов
    figure;
    plot(z, I);
    xlabel('z (м)');
    ylabel('Интенсивность');
    title(['Распределение интенсивности для NA = ',num2str(NA_val)]);

end
