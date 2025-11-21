q = 1.6e-19;
m = 9.11e-31;
B = 0.000000001;
E = 1;

k1 = q * B / m;
k2 = q * E / m;

omega_c = q * B / m;  // Циклотронная частота
a_E = q * E / m;  // Ускорение в электрическом поле

// Магнитное поле
function dYdt = fB(t, Y)
    x = Y(1);
    y = Y(2);
    vx = Y(3);
    vy = Y(4);

    dxdt = vx;
    dydt = vy;
    dvxdt = omega_c * vy;
    dvydt = -omega_c * vx;

    dYdt = [dxdt; dydt; dvxdt; dvydt];
endfunction

// Электро поле
function dYdt = fE(t, Y)
    x = Y(1);
    y = Y(2);
    vx = Y(3);
    vy = Y(4);

    dxdt = vx;
    dydt = vy;
    dvxdt = 0;
    dvydt = -a_E;

    dYdt = [dxdt; dydt; dvxdt; dvydt];
endfunction

// Электро-магнитное
function dYdt = fEB(t, Y)
    x = Y(1);
    y = Y(2);
    vx = Y(3);
    vy = Y(4);

    dxdt = vx;
    dydt = vy;
    dvxdt = omega_c * vy + 0;
    dvydt = -omega_c * vx + a_E;

    dYdt = [dxdt; dydt; dvxdt; dvydt];
endfunction

Y0 = [0; 0; 1; 0]; // [x0; y0; vx0; vy0]
t0 = 0;
t = 0:0.01:30;

YB = ode(Y0, t0, t, fB);
YE = ode(Y0, t0, t, fE);
YEB = ode(Y0, t0, t, fEB);

scf();
subplot(1,3,1);
plot(YB(1,:), YB(2,:), 'b-');
title('Магнитное поле');
xlabel('x'); ylabel('y');
xgrid;

subplot(1,3,2);
plot(YE(1,:), YE(2,:), 'r-');
title('Электрическое поле');
xlabel('x'); ylabel('y');
xgrid;

subplot(1,3,3);
plot(YEB(1,:), YEB(2,:), 'black-');
title('Электромагнитное поле');
xlabel('x'); ylabel('y');
xgrid;
