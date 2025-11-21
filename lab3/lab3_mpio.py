import numpy as np
import matplotlib.pyplot as plt
import os
import math

# Создаём папку для результатов
os.makedirs("results", exist_ok=True)

def rk4_step(func, y, t, dt):
    """Один шаг метода Рунге–Кутты 4-го порядка."""
    k1 = np.array(func(t, y))
    k2 = np.array(func(t + dt/2, y + dt/2 * k1))
    k3 = np.array(func(t + dt/2, y + dt/2 * k2))
    k4 = np.array(func(t + dt, y + dt * k3))
    return y + dt / 6 * (k1 + 2*k2 + 2*k3 + k4)

def solve_ode(func, y0, t_span, dt):
    """Решает систему ОДУ методом РК4 на интервале t_span с шагом dt."""
    t0, tf = t_span
    N = int((tf - t0) / dt) + 1
    t_vals = np.linspace(t0, tf, N)
    y_vals = [np.array(y0, dtype=float)]
    y = np.array(y0, dtype=float)
    for i in range(1, N):
        y = rk4_step(func, y, t_vals[i-1], dt)
        y_vals.append(y.copy())
    return np.array(t_vals), np.array(y_vals)

# ----------------------------
# ЗАДАНИЕ 1: Осцилляторы
# ----------------------------

def linear_oscillator(t, y, omega=1.0, gamma=0.0):
    x, v = y
    dxdt = v
    dvdt = -omega**2 * x - gamma * v
    return [dxdt, dvdt]

def nonlinear_oscillator(t, y, omega=1.0, gamma=0.0):
    x, v = y
    dxdt = v
    dvdt = -omega**2 * np.sin(x) - gamma * v
    return [dxdt, dvdt]

def task1():
    omega = 1.0
    t_span = (0, 30)
    dt = 0.01
    init_conditions = [
        [1.0, 0.0],
        [2.0, 0.5],
        [0.5, -1.0]
    ]
    gammas = [0.0, 0.4]

    plt.figure(figsize=(12, 5))
    for i, gamma in enumerate(gammas):
        plt.subplot(1, 2, i+1)
        for y0 in init_conditions:
            sol = lambda t, y: linear_oscillator(t, y, omega, gamma)
            t, y = solve_ode(sol, y0, t_span, dt)
            plt.plot(y[:, 0], y[:, 1], label=f'x0={y0[0]}, v0={y0[1]}')
        plt.title(f'Линейный осциллятор (γ = {gamma})')
        plt.xlabel('x')
        plt.ylabel('v')
        plt.grid(True)
        plt.legend()
    plt.tight_layout()
    plt.savefig("results/task1_linear.png", dpi=150)
    plt.show()

    plt.figure(figsize=(12, 5))
    for i, gamma in enumerate(gammas):
        plt.subplot(1, 2, i+1)
        for y0 in init_conditions:
            sol = lambda t, y: nonlinear_oscillator(t, y, omega, gamma)
            t, y = solve_ode(sol, y0, t_span, dt)
            plt.plot(y[:, 0], y[:, 1], label=f'x0={y0[0]}, v0={y0[1]}')
        plt.title(f'Нелинейный осциллятор (γ = {gamma})')
        plt.xlabel('x')
        plt.ylabel('v')
        plt.grid(True)
        plt.legend()
    plt.tight_layout()
    plt.savefig("results/task1_nonlinear.png", dpi=150)
    plt.show()

# ----------------------------
# ЗАДАНИЕ 2: Брюсселятор
# ----------------------------

def brusselator(t, y, mu=1.0):
    u1, u2 = y
    du1dt = 1 - (mu + 1) * u1 + u1**2 * u2
    du2dt = mu * u1 - u1**2 * u2
    return [du1dt, du2dt]

def task2():
    t_span = (0, 50)
    dt = 0.01
    init_conditions = [
        [0.5, 0.5],
        [1.5, 2.0],
        [2.0, 1.0]
    ]
    mus = [0.5, 2.5]

    plt.figure(figsize=(12, 5))
    for i, mu in enumerate(mus):
        plt.subplot(1, 2, i+1)
        for y0 in init_conditions:
            sol = lambda t, y: brusselator(t, y, mu)
            t, y = solve_ode(sol, y0, t_span, dt)
            plt.plot(y[:, 0], y[:, 1], label=f'u1(0)={y0[0]}, u2(0)={y0[1]}')
        plt.title(f'Брюсселятор (μ = {mu})')
        plt.xlabel('u1')
        plt.ylabel('u2')
        plt.grid(True)
        plt.legend()
    plt.tight_layout()
    plt.savefig("results/task2_brusselator.png", dpi=150)
    plt.show()

# ----------------------------
# ЗАДАНИЕ 3: Модель хищник-жертва
# ----------------------------

def lotka_volterra(t, y, A1=0.1, D2=1.0, D12=0.1, A21=0.1, B=0.0):
    x, y_ = y
    dxdt = A1 * x - D12 * x * y_ - B * x**2 + math.sin(t/3) ** 2 * 3
    dydt = -D2 * y_ + A21 * x * y_ + math.sin(t/3) ** 2 * 3
    print(dxdt)
    return [dxdt, dydt]

def task3():
    t_span = (0, 100)
    dt = 0.01
    init_conditions = [
        [10, 5],
        [20, 10],
        [5, 15],
    ]
    Bs = [0.0, 0.001]  # B=0 — классическая модель; B>0 — с конкуренцией

    plt.figure(figsize=(12, 5))
    for i, B in enumerate(Bs):
        plt.subplot(1, 2, i+1)
        for y0 in init_conditions:
            sol = lambda t, y: lotka_volterra(t, y, B=B)
            t, y = solve_ode(sol, y0, t_span, dt)
            plt.plot(y[:, 0], y[:, 1], label=f'x0={y0[0]}, y0={y0[1]}')
        title = 'Хищник-жертва (без конкуренции)' if B == 0 else f'Хищник-жертва (B = {B})'
        plt.title(title)
        plt.xlabel('Жертвы (x)')
        plt.ylabel('Хищники (y)')
        plt.grid(True)
        plt.legend()
    plt.tight_layout()
    plt.savefig("results/task3_predator_prey.png", dpi=150)
    plt.show()

# ----------------------------
# Запуск всех заданий
# ----------------------------

if __name__ == "__main__":
    print("Выполняется задание 1: Осцилляторы...")
    # task1()
    print("Выполняется задание 2: Брюсселятор...")
    # task2()
    print("Выполняется задание 3: Модель хищник-жертва...")
    task3()
    print("Готово!")
