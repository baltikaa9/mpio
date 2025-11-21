K = 4 * %pi^2;
M = 1;
m1 = 0.001;
m2 = 0.01;

function dYdt = ThreeBodyRHP(t, Y)
    x1 = Y(1);
    y1 = Y(2);
    vx1 = Y(3);
    vy1 = Y(4);
    x2 = Y(5);
    y2 = Y(6);
    vx2 = Y(7);
    vy2 = Y(8);

    r1 = sqrt(x1^2 + y1^2);
    r2 = sqrt(x2^2 + y2^2);
    dx21 = x2 - x1;
    dy21 = y2 - y1;
    r21 = sqrt(dx21^2 + dy21^2);

    dx1dt = vx1;
    dy1dt = vy1;
    dx2dt = vx2;
    dy2dt = vy2;

    dvx1dt = -K * x1 / r1^3 + K * m2 * dx21 / r21^3;
    dvy1dt = -K * y1 / r1^3 + K * m2 * dy21 / r21^3;

    dvx2dt = -K * x2 / r2^3 - K * m1 * dx21 / r21^3;
    dvy2dt = -K * y2 / r2^3 - K * m1 * dy21 / r21^3;

    dYdt = [dx1dt; dy1dt; dvx1dt; dvy1dt; dx2dt; dy2dt; dvx2dt; dvy2dt];
endfunction

function main(x1_0, y1_0, vx1_0, vy1_0, x2_0, y2_0, vx2_0, vy2_0, tmax)
    disp('=== Задача трех тел ===')
    disp('Масса звезды M = ' + string(M) + ' M☉')
    disp('Масса планеты 1 m1 = ' + string(m1) + ' M☉')
    disp('Масса планеты 2 m2 = ' + string(m2) + ' M☉')
    disp(' ')
    disp('Начальные условия планеты 1:')
    disp('x1 = ' + string(x1_0) + ' а.е., y1 = ' + string(y1_0) + ' а.е.')
    disp('vx1 = ' + string(vx1_0) + ', vy1 = ' + string(vy1_0))
    disp(' ')
    disp('Начальные условия планеты 2:')
    disp('x2 = ' + string(x2_0) + ' а.е., y2 = ' + string(y2_0) + ' а.е.')
    disp('vx2 = ' + string(vx2_0) + ', vy2 = ' + string(vy2_0))

    Y0 = [x1_0; y1_0; vx1_0; vy1_0; x2_0; y2_0; vx2_0; vy2_0];
    t0 = 0;
    t = 0:0.01:tmax;

    Y = ode(Y0, t0, t, ThreeBodyRHP);

    x1_array = Y(1, :);
    y1_array = Y(2, :);
    vx1_array = Y(3, :);
    vy1_array = Y(4, :);
    x2_array = Y(5, :);
    y2_array = Y(6, :);
    vx2_array = Y(7, :);
    vy2_array = Y(8, :);

    E = zeros(1, length(t));
    L = zeros(1, length(t));

    for i = 1:length(t)
        x1 = x1_array(i);
        y1 = y1_array(i);
        vx1 = vx1_array(i);
        vy1 = vy1_array(i);
        x2 = x2_array(i);
        y2 = y2_array(i);
        vx2 = vx2_array(i);
        vy2 = vy2_array(i);

        r1 = sqrt(x1^2 + y1^2);
        r2 = sqrt(x2^2 + y2^2);
        r21 = sqrt((x2-x1)^2 + (y2-y1)^2);

        E_kin = 0.5 * m1 * (vx1^2 + vy1^2) + 0.5 * m2 * (vx2^2 + vy2^2);
        E_pot = -K * m1 / r1 - K * m2 / r2 - K * m1 * m2 / r21;
        E(i) = E_kin + E_pot;

        L1 = m1 * (x1 * vy1 - y1 * vx1);
        L2 = m2 * (x2 * vy2 - y2 * vx2);
        L(i) = L1 + L2;
    end

    E0 = E(1);
    L0 = L(1);

    disp(' ')
    disp('=== Сохранение энергии ===')
    disp('Начальная энергия E0 = ' + string(E0))
    disp('Конечная энергия E = ' + string(E($)))
    disp('Относительное изменение = ' + string(abs(E($) - E0) / abs(E0) * 100) + ' %')

    disp(' ')
    disp('=== Сохранение момента импульса ===')
    disp('Начальный момент импульса L0 = ' + string(L0))
    disp('Конечный момент импульса L = ' + string(L($)))
    disp('Относительное изменение = ' + string(abs(L($) - L0) / abs(L0) * 100) + ' %')

    // Формируем строку с НУ
    ic_str = 'НУ: r1=(' + string(x1_0) + ',' + string(y1_0) + ') r2=(' + string(x2_0) + ',' + string(y2_0) + ')';
    
    // Построение графиков
    fig = scf();
    fig.figure_name = ic_str;  // Заголовок окна с НУ

    subplot(1,3,1);
    plot(x1_array, y1_array, 'b-', 'LineWidth', 1);
    plot(x2_array, y2_array, 'r-', 'LineWidth', 1);
    plot(0, 0, 'o', 'MarkerSize', 15, 'MarkerFaceColor', '#FFD700', 'MarkerEdgeColor', '#FFD700');
    plot(x1_0, y1_0, 'bs', 'MarkerSize', 8);
    plot(x2_0, y2_0, 'rs', 'MarkerSize', 8);
    legend(['Планета 1'; 'Планета 2'; 'Звезда'], 'in_upper_left');
    title('Траектории планет');
    xlabel('x (а.е.)'); ylabel('y (а.е.)');
    xgrid;

    subplot(1,3,2);
    plot(t, E, 'b-');
    title('Полная энергия');
    xlabel('Время (годы)'); ylabel('Энергия E');
    xgrid;

    subplot(1,3,3);
    plot(t, L, 'r-');
    title('Момент импульса');
    xlabel('Время (годы)'); ylabel('Момент импульса L');
    xgrid;
endfunction

disp(' ')
disp('========================================')
disp('Пример 1: Стабильная конфигурация')
disp('========================================')

x1_0 = 1;
y1_0 = 0;
vx1_0 = 0;
vy1_0 = sqrt(K / x1_0);

x2_0 = 2.5;
y2_0 = 0;
vx2_0 = 0;
vy2_0 = sqrt(K / x2_0);

main(x1_0, y1_0, vx1_0, vy1_0, x2_0, y2_0, vx2_0, vy2_0, 100);

disp('========================================')
disp('Пример 2')
disp('========================================')

x1_0 = 1;
y1_0 = 0;
vx1_0 = 0;
vy1_0 = sqrt(K / x1_0);

x2_0 = 1.5;
y2_0 = 0;
vx2_0 = 0;
vy2_0 = sqrt(K / x2_0);

main(x1_0, y1_0, vx1_0, vy1_0, x2_0, y2_0, vx2_0, vy2_0, 20);

disp(' ')
disp('========================================')
disp('Пример 3: Взаимодействующие планеты')
disp('========================================')

x1_0 = 1;
y1_0 = 0;
vx1_0 = 0;
vy1_0 = sqrt(K / x1_0);

x2_0 = 1.3;
y2_0 = 0;
vx2_0 = 0;
vy2_0 = sqrt(K / x2_0);

main(x1_0, y1_0, vx1_0, vy1_0, x2_0, y2_0, vx2_0, vy2_0, 50);

disp(' ')
disp('========================================')
disp('Пример 4')
disp('========================================')

x1_0 = 1;
y1_0 = 0;
vx1_0 = 0;
vy1_0 = sqrt(K / x1_0);

x2_0 = -1.2;
y2_0 = 0;
vx2_0 = 0;
vy2_0 = sqrt(K / abs(x2_0));

main(x1_0, y1_0, vx1_0, vy1_0, x2_0, y2_0, vx2_0, vy2_0, 50);

disp(' ')
disp('========================================')
disp('Пример 5: Столкновение')
disp('========================================')

x1_0 = 1;
y1_0 = 0;
vx1_0 = 0;
vy1_0 = sqrt(K / x1_0);

x2_0 = -1;
y2_0 = 0;
vx2_0 = 0;
vy2_0 = sqrt(K / abs(x2_0));

main(x1_0, y1_0, vx1_0, vy1_0, x2_0, y2_0, vx2_0, vy2_0, 10);
