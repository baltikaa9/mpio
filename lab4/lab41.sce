K = 4 * %pi^2;

function dYdt = RHP(t, Y)
    x = Y(1);
    y = Y(2);
    vx = Y(3);
    vy = Y(4);

    dxdt = vx;
    dydt = vy;
    dvxdt = -K * x / (x^2 + y^2)^1.5;
    dvydt = -K * y / (x^2 + y^2)^1.5;

    dYdt = [dxdt; dydt; dvxdt; dvydt];
endfunction

function main(x0, y0, vx0, vy0, tmax, filename)
    disp('Начальные условия:')
    disp('x0 = ' + string(x0))
    disp('y0 = ' + string(y0))
    disp('vx0 = ' + string(vx0))
    disp('vy0 = ' + string(vy0))

    Y0 = [x0; y0; vx0; vy0];
    t0 = 0;
    t = 0:0.001:tmax;

    Y = ode(Y0, t0, t, RHP);

    // Вычисление эксцентриситета орбиты
    // Используем начальные условия: r = (x0, y0), v = (vx0, vy0)
    r0_vec = [x0; y0; 0];
    v0_vec = [vx0; vy0; 0];
    r0 = sqrt(x0^2 + y0^2);
    v0 = sqrt(vx0^2 + vy0^2);

    // Вектор момента импульса: L = r × v
    L_vec = cross(r0_vec, v0_vec);
    L = sqrt(L_vec(1)^2 + L_vec(2)^2 + L_vec(3)^2);

    // Полная энергия в начальный момент
    E_orbit = v0^2 / 2 - K / r0;

    // Большая полуось: a = -K / (2*E)
    a_orbit = -K / (2 * E_orbit);

    // Эксцентриситет: e = sqrt(1 + 2*E*L^2/K^2)
    eccentricity = sqrt(1 + 2 * E_orbit * L^2 / K^2);

    // Перицентр и апоцентр
    r_peri = a_orbit * (1 - eccentricity);
    r_apo = a_orbit * (1 + eccentricity);

    // Вычисление полной энергии и момента импульса на каждом шаге
    E = zeros(1, length(t));
    L_array = zeros(1, length(t));
    for i = 1:length(t)
        x = Y(1, i);
        y = Y(2, i);
        vx = Y(3, i);
        vy = Y(4, i);

        // Кинетическая энергия: E_k = (vx^2 + vy^2) / 2
        E_kinetic = (vx^2 + vy^2) / 2;

        // Потенциальная энергия: E_p = -K / r = -K / sqrt(x^2 + y^2)
        r = sqrt(x^2 + y^2);
        E_potential = -K / r;

        // Полная энергия
        E(i) = E_kinetic + E_potential;

        // Момент импульса: L = |r × v| = |x*vy - y*vx|
        L_array(i) = abs(x * vy - y * vx);
    end

    // Начальная энергия (для сравнения)
    E0 = E(1);
    disp(' ')
    disp('=== Энергия ===')
    disp('Начальная полная энергия E0 = ' + string(E0))
    disp('Конечная полная энергия E = ' + string(E($)))
    disp('Изменение энергии dE = ' + string(E($) - E0))
    disp('Относительное изменение = ' + string(abs(E($) - E0) / abs(E0) * 100) + ' %')

    // Момент импульса (для сравнения)
    L0 = L_array(1);
    disp(' ')
    disp('=== Момент импульса ===')
    disp('Начальный момент импульса L0 = ' + string(L0))
    disp('Конечный момент импульса L = ' + string(L_array($)))
    disp('Изменение момента импульса dL = ' + string(L_array($) - L0))
    disp('Относительное изменение = ' + string(abs(L_array($) - L0) / abs(L0) * 100) + ' %')

    disp(' ')
    disp('=== Параметры орбиты ===')
    disp('Эксцентриситет e = ' + string(eccentricity))
    if eccentricity < 0.01 then
        disp('Тип орбиты: круговая')
    elseif eccentricity < 1 then
        disp('Тип орбиты: эллиптическая')
    elseif abs(eccentricity - 1) < 0.01 then
        disp('Тип орбиты: параболическая')
    else
        disp('Тип орбиты: гиперболическая')
    end
    disp('Большая полуось a = ' + string(a_orbit) + ' а.е.')
    disp('Перицентр r_p = ' + string(r_peri) + ' а.е.')
    disp('Апоцентр r_a = ' + string(r_apo) + ' а.е.')
    disp('Момент импульса L = ' + string(L))

    // Находим момент, когда планета пересекает ось X второй раз (y ≈ 0, x > 0)
    T_numerical = 0;
    for i = 2:length(t)
        if Y(2, i-1) < 0 && Y(2, i) >= 0 && Y(1, i) > 0 then
            // Линейная интерполяция для точного определения момента пересечения
            T_numerical = t(i-1) + (0 - Y(2, i-1)) * (t(i) - t(i-1)) / (Y(2, i) - Y(2, i-1));
            break;
        end
    end

    // Теоретический период по 3-му закону Кеплера (используем вычисленную большую полуось)
    T_theoretical = 2 * %pi * sqrt(a_orbit^3 / K);

    disp(' ')
    disp('=== Период обращения ===')
    disp('Численный период T = ' + string(T_numerical) + ' лет')
    disp('Теоретический период (3 закон Кеплера) T = ' + string(T_theoretical) + ' лет')
    disp('Относительная погрешность = ' + string(abs(T_numerical - T_theoretical) / T_theoretical * 100) + ' %')
    disp(' ')

    // Формируем строку с НУ
    ic_str = 'НУ: r=(' + string(x0) + ',' + string(y0) + '), v=(' + string(vx0) + ',' + string(vy0) + ')';

    fig = scf();
    fig.figure_name = ic_str;

    subplot(1,3,1);
    plot(Y(1,:), Y(2,:), 'black-');
    plot(0, 0, 'o', 'MarkerSize', 15, 'MarkerFaceColor', '#FFD700', 'MarkerEdgeColor', '#FFD700');
    title('Движение по орбите');
    xlabel('x (а.е.)'); ylabel('y (а.е.)');
    xgrid;

    subplot(1,3,2);
    plot(t, E, 'b-');
    title('Полная энергия системы от времени');
    xlabel('Время (годы)'); ylabel('Полная энергия E');
    xgrid;

    subplot(1,3,3);
    plot(t, L_array, 'r-');
    title('Момент импульса от времени');
    xlabel('Время (годы)'); ylabel('Момент импульса L');
    xgrid;

    // Сохранение графика
    xs2png(gcf(), 'img/' + filename);
    disp('График сохранен в файл img/' + filename);
endfunction

// Круговая орбита (e = 0)
x0 = 1;
y0 = 0;
vx0 = 0;
vy0 = sqrt(K / x0);  // Первая космическая скорость

main(x0, y0, vx0, vy0, 30, 'orbit_circular.png');

// Эллиптическая орбита

// Вариант 1: Слабо вытянутый эллипс (e ≈ 0.3)
x0 = 1;
y0 = 0;
vx0 = 0;
vy0 = 1.2 * sqrt(K / x0);  // На 20% больше круговой

main(x0, y0, vx0, vy0, 30, 'orbit_elliptic_weak.png')

// Вариант 2: Умеренный эллипс (e ≈ 0.5)
x0 = 1;
y0 = 0;
vx0 = 0;
vy0 = 1.4 * sqrt(K / x0);  // На 40% больше круговой

main(x0, y0, vx0, vy0, 200, 'orbit_elliptic_moderate.png');

// Параболическая орбита (e = 1)
x0 = 1;
y0 = 0;
vx0 = 0;
vy0 = sqrt(2 * K / x0);  // Вторая космическая скорость

main(x0, y0, vx0, vy0, 100, 'orbit_parabolic.png');

// Гиперболическая орбита (e > 1)
x0 = 1;
y0 = 0;
vx0 = 0;
vy0 = 1.02 * sqrt(2 * K / x0);

main(x0, y0, vx0, vy0, 50, 'orbit_hyperbolic_close.png');

// Гиперболическая орбита - дальний подлет
x0 = 5;
y0 = -20;
vx0 = 0;
vy0 = 1.05 * sqrt(2 * K / x0);

main(x0, y0, vx0, vy0, 50, 'orbit_hyperbolic_distant.png');
