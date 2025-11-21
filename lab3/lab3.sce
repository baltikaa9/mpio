function dYdt = linear_oscillator(t, Y)
    x = Y(1);
    y = Y(2);

    dxdt = y;
    dydt = -omega^2 * x + gamma * y;

    dYdt = [dxdt; dydt];
endfunction

function dYdt = nonlinear_oscillator(t, Y)
    x = Y(1);
    y = Y(2);

    dxdt = y;
    dydt = -(omega^2) * sin(x) - gamma * y;

    dYdt = [dxdt; dydt];
endfunction

function dYdt = brusselator(t, Y)
    u1 = Y(1);
    u2 = Y(2);

    du1dt = -(gamma + 1) * u1 + u1^2 * u2 + 1;  // вместо mu - gamma
    du2dt = gamma * u1 - u1^2 * u2;

    dYdt = [du1dt; du2dt];
endfunction

function dYdt = lotka_volterra(t, Y)
    x = Y(1);
    y = Y(2);

    dxdt = A1 * x - D1 * x * y - B * x^2 + sin(t/3)^2 * 3;
    dydt = -D2 * y + A2 * x * y + sin(t/3)^2 * 3;

    dYdt = [dxdt; dydt];
endfunction

Y0s = [[2; 0], [0; 2], [2; 2]]; // [x0; y0]
t0 = 0;
t = 0:0.01:10;

omega = 1;
gammas = [0, 0.3, 0.5];

colors = ['r', 'black', 'b'];

// Линейный осциллятор
scf();
for i = 1:length(gammas)
    gamma = gammas(i);
    subplot(1,length(gammas),i);
    labels = [];

    for j = 1:size(Y0s, 2)
        Y = ode(Y0s(:,j), t0, t, linear_oscillator);

        plot(Y(1,:), Y(2,:), colors(j));
        labels(j) = msprintf('x = %g, v = %g', Y0s(1,j), Y0s(2,j));
    end
    title(msprintf('Нелинейный осциллятор, $\gamma$ = %g', gamma));
    xlabel('x'); ylabel('v');
    legend(labels);
    xgrid;
end
xs2png(gcf(), 'img/linear_oscillator.png');
disp('График линейного осциллятора сохранен в linear_oscillator.png');

// Нелинейный осциллятор

//Y0s = [[1.0; 0.0], [2.0; 0.5], [0.5; -1.0]];

scf();
for i = 1:length(gammas)
    gamma = gammas(i);
    subplot(1,length(gammas),i);
    labels = [];

    for j = 1:size(Y0s, 2)
        Y = ode(Y0s(:,j), t0, t, nonlinear_oscillator);

        plot(Y(1,:), Y(2,:), colors(j));
        labels(j) = msprintf('x = %g, v = %g', Y0s(1,j), Y0s(2,j));
    end
    title(msprintf('Нелинейный осциллятор, $\gamma$ = %g', gamma));
    xlabel('x'); ylabel('v');
    legend(labels);
    xgrid;
end
xs2png(gcf(), 'img/nonlinear_oscillator.png');
disp('График нелинейного осциллятора сохранен в nonlinear_oscillator.png');

// Брюсселятор

//Y0s = [[1.0; 0.0], [2.0; 0.5], [0.5; -1.0]];
gammas = [0, 1, 2];

scf();
for i = 1:length(gammas)
    gamma = gammas(i);  // вместо mu
    subplot(1,length(gammas),i);
    labels = [];

    for j = 1:size(Y0s, 2)
        Y = ode(Y0s(:,j), t0, t, brusselator);

        plot(Y(1,:), Y(2,:), colors(j) + '.', 'MarkerSize', 3);
        labels(j) = msprintf('u1(0) = %g, u2(0) = %g', Y0s(1,j), Y0s(2,j));
    end
    title(msprintf('Брюсселятор, $\mu$ = %g', gamma));
    xlabel('u1'); ylabel('u2');
    legend(labels);
    xgrid;
end
xs2png(gcf(), 'img/brusselator.png');
disp('График брюсселятора сохранен в brusselator.png');

// Лотка-Вольтерра

Y0s = [[10; 5],[10; 10],[10; 20],[10; 30]]; // [x0; y0]
A1 = 0.1;
A2 = 0.1;
D1 = 0.1;
D2 = 1;

Bs = [0, 0.5, 1];

t = 0:0.01:100;

scf();
for i = 1:length(Bs)
    B = Bs(i);
    subplot(1,length(Bs),i);
    labels = [];

    for j = 1:size(Y0s, 2)
        Y = ode(Y0s(:,j), t0, t, lotka_volterra);

        plot(Y(1,:), Y(2,:), 'black.', 'MarkerSize', 3);
        labels(j) = msprintf('x = %g, y = %g', Y0s(1,j), Y0s(2,j));
    end
    title(msprintf('Лотка-Вольтерра, B = %g', B));
    xlabel('x'); ylabel('y');
    legend(labels);
    xgrid;
end
xs2png(gcf(), 'img/lotka_volterra.png');
disp('График Лотка-Вольтерра сохранен в lotka_volterra.png');
