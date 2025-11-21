import matplotlib.pyplot as plt
import numpy as np

# --- ЧАСТЬ 1: ОСНОВНЫЕ ФУНКЦИИ И ИНТЕГРАТОР ---

# Константа K в безразмерной системе единиц (а.е., годы)
K = 4 * np.pi ** 2


def rk4_step(func, y, t, dt, **kwargs):
    """
    Один шаг интегрирования методом Рунге-Кутты 4-го порядка.

    Аргументы:
        func: функция, вычисляющая производные dy/dt.
        y: текущее состояние системы (вектор).
        t: текущее время.
        dt: шаг по времени.
        **kwargs: дополнительные аргументы для функции func (например, K).

    Возвращает:
        y_next: состояние системы в момент времени t + dt.
    """
    k1 = dt * func(t, y, **kwargs)
    k2 = dt * func(t + 0.5 * dt, y + 0.5 * k1, **kwargs)
    k3 = dt * func(t + 0.5 * dt, y + 0.5 * k2, **kwargs)
    k4 = dt * func(t + dt, y + k3, **kwargs)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def two_body_derivatives(t, y, K_val):
    """
    Вычисляет производные для задачи двух тел (Солнце + 1 планета).
    Вектор состояния y = [x, y, vx, vy].
    """
    x, y_pos, vx, vy = y
    r_cubed = (x ** 2 + y_pos ** 2) ** 1.5

    # dx/dt = vx
    # dy/dt = vy
    # dvx/dt = ax
    # dvy/dt = ay

    ax = -K_val * x / r_cubed
    ay = -K_val * y_pos / r_cubed

    return np.array([vx, vy, ax, ay])


def simulate(func, y0, t_span, dt, **kwargs):
    """
    Запускает симуляцию для решения системы ОДУ.
    """
    t_values = np.arange(t_span[0], t_span[1], dt)
    y_values = np.zeros((len(t_values), len(y0)))
    y_values[0] = y0

    for i in range(1, len(t_values)):
        y_values[i] = rk4_step(func, y_values[i - 1], t_values[i - 1], dt, **kwargs)

    return t_values, y_values


# --- ЧАСТЬ 2: РЕШЕНИЕ ЗАДАНИЙ ---

def plot_orbit_and_energy(t, results, title, K_val):
    """Вспомогательная функция для отрисовки графиков."""
    x = results[:, 0]
    y = results[:, 1]
    vx = results[:, 2]
    vy = results[:, 3]

    # Вычисление энергии и момента импульса
    r = np.sqrt(x ** 2 + y ** 2)
    v_sq = vx ** 2 + vy ** 2
    # Масса планеты m=1 в безразмерных единицах
    energy = 0.5 * v_sq - K_val / r
    angular_momentum = x * vy - y * vx

    fig, axs = plt.subplots(1, 3, figsize=(21, 6))
    fig.suptitle(title, fontsize=16)

    # 1. График орбиты
    axs[0].plot(x, y, label='Траектория планеты')
    axs[0].plot(0, 0, 'yo', markersize=10, label='Солнце')
    axs[0].set_title('Орбита')
    axs[0].set_xlabel('x (а.е.)')
    axs[0].set_ylabel('y (а.е.)')
    axs[0].axis('equal')
    axs[0].grid(True)
    axs[0].legend(loc='upper right')

    # 2. График сохранения энергии
    axs[1].plot(t, (energy - energy[0]) / energy[0])
    axs[1].set_title('Сохранение энергии')
    axs[1].set_xlabel('Время (годы)')
    axs[1].set_ylabel('Относительное изменение E: (E(t)-E₀)/E₀')
    axs[1].grid(True)

    # 3. График сохранения момента импульса
    axs[2].plot(t, (angular_momentum - angular_momentum[0]) / angular_momentum[0])
    axs[2].set_title('Сохранение момента импульса')
    axs[2].set_xlabel('Время (годы)')
    axs[2].set_ylabel('Относительное изменение L: (L(t)-L₀)/L₀')
    axs[2].grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

    # Проверка третьего закона Кеплера
    # Найдем период, отследив полный оборот
    y_cross_indices = np.where((y[:-1] < 0) & (y[1:] >= 0))[0]
    if len(y_cross_indices) > 1:
        period = t[y_cross_indices[1]] - t[y_cross_indices[0]]
        a = (np.max(r) + np.min(r)) / 2  # Большая полуось
        kepler_const = period ** 2 / a ** 3
        print(f"Большая полуось a = {a:.4f} а.е.")
        print(f"Вычисленный период T = {period:.4f} лет")
        print(f"Проверка 3-го закона Кеплера: T²/a³ = {kepler_const:.4f}")
        print(f"Теоретическое значение: 4π²/K = {4 * np.pi ** 2 / K_val:.4f}")


### ЗАДАНИЕ 5.2.1: Проверка программы для круговых орбит
print("--- Задание 5.2.1: Круговая орбита ---")
# Начальные условия для Земли (круговая орбита)
# x₀=1, y₀=0, vx₀=0. Скорость voy для круговой орбиты: v = sqrt(K/r)
x0_circ = 1.0
v0y_circ = np.sqrt(K / x0_circ)
y0_circ = np.array([x0_circ, 0.0, 0.0, v0y_circ])

# Симуляция на 5 лет с шагом 0.001 года
t_circ, results_circ = simulate(two_body_derivatives, y0_circ, (0, 5), 0.001, K_val=K)
plot_orbit_and_energy(t_circ, results_circ, 'Модель круговой орбиты (Земля)', K)

### ЗАДАНИЕ 5.2.2: Проверка программы для эллиптических орбит
print("\n--- Задание 5.2.2: Эллиптическая орбита ---")
# Начальные условия для эллиптической орбиты. Возьмем скорость меньше,
# чем для круговой, чтобы получить эллипс.
x0_ell = 1.0  # Начинаем в афелии (самой дальней точке)
v0y_ell = 0.3 * np.sqrt(K / x0_ell)  # Скорость на 20% меньше
y0_ell = np.array([x0_ell, 0.0, 0.0, v0y_ell])
print(y0_ell)

# Симуляция на 5 лет
t_ell, results_ell = simulate(two_body_derivatives, y0_ell, (0, 5), 0.001, K_val=K)
plot_orbit_and_energy(t_ell, results_ell, 'Модель эллиптической орбиты', K)

### ЗАДАНИЕ 5.2.3: Маленькая Солнечная система (Задача трех тел)
print("\n--- Задание 5.2.3: Маленькая Солнечная система ---")


def three_body_derivatives(t, y, K_val, m1_ratio, m2_ratio):
    """
    Вычисляет производные для задачи трех тел (Солнце + 2 планеты).
    Вектор состояния y = [x1, y1, vx1, vy1, x2, y2, vx2, vy2].
    """
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = y

    # Расстояния
    r1_sq = x1 ** 2 + y1 ** 2
    r2_sq = x2 ** 2 + y2 ** 2
    r12_sq = (x1 - x2) ** 2 + (y1 - y2) ** 2

    r1_cubed = r1_sq ** 1.5
    r2_cubed = r2_sq ** 1.5
    r12_cubed = r12_sq ** 1.5

    # Ускорения для планеты 1 (от Солнца и от планеты 2)
    ax1 = -K_val * x1 / r1_cubed - K_val * m2_ratio * (x1 - x2) / r12_cubed
    ay1 = -K_val * y1 / r1_cubed - K_val * m2_ratio * (y1 - y2) / r12_cubed

    # Ускорения для планеты 2 (от Солнца и от планеты 1)
    ax2 = -K_val * x2 / r2_cubed - K_val * m1_ratio * (x2 - x1) / r12_cubed
    ay2 = -K_val * y2 / r2_cubed - K_val * m1_ratio * (y2 - y1) / r12_cubed

    return np.array([vx1, vy1, ax1, ay1, vx2, vy2, ax2, ay2])


# Параметры из методички
m1_ratio = 0.001  # Отношение массы планеты 1 к массе Солнца
m2_ratio = 0.01  # Отношение массы планеты 2 к массе Солнца

# Начальные условия
x10 = 1.0
y10 = 0.0
vx10 = 0.0
vy10 = np.sqrt(K / x10)  # Начальная скорость для круговой орбиты

x20 = 4.0 ** (1 / 3)
y20 = 0.0
vx20 = 0.0
vy20 = np.sqrt(K / x20)  # Начальная скорость для круговой орбиты

y0_3body = np.array([x10, y10, vx10, vy10, x20, y20, vx20, vy20])

# Симуляция на 20 лет
t_3body, results_3body = simulate(three_body_derivatives, y0_3body, (0, 20), 0.001, K_val=K, m1_ratio=m1_ratio,
                                  m2_ratio=m2_ratio)

# --- Анализ результатов для задачи трех тел ---
x1, y1, vx1, vy1, x2, y2, vx2, vy2 = results_3body.T

# Энергии планет по отдельности (без учета взаимодействия)
r1 = np.sqrt(x1 ** 2 + y1 ** 2)
r2 = np.sqrt(x2 ** 2 + y2 ** 2)
energy1 = 0.5 * m1_ratio * (vx1 ** 2 + vy1 ** 2) - K * m1_ratio / r1
energy2 = 0.5 * m2_ratio * (vx2 ** 2 + vy2 ** 2) - K * m2_ratio / r2

# Полная энергия системы (включая потенциальную энергию взаимодействия планет)
r12 = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
pe_interaction = -K * m1_ratio * m2_ratio / r12
total_energy = energy1 + energy2 + pe_interaction

# Построение графиков для задачи трех тел
fig, axs = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle('Задание 5.2.3: Модель "Маленькой Солнечной системы"', fontsize=16)

# 1. График орбит
axs[0].plot(x1, y1, label='Планета 1 (m/M=0.001)')
axs[0].plot(x2, y2, label=f'Планета 2 (m/M={m2_ratio})')
axs[0].plot(0, 0, 'yo', markersize=10, label='Солнце')
axs[0].set_title('Орбиты планет')
axs[0].set_xlabel('x (а.е.)')
axs[0].set_ylabel('y (а.е.)')
axs[0].axis('equal')
axs[0].grid(True)
axs[0].legend(loc='upper right')

# 2. График сохранения энергии
axs[1].plot(t_3body, (energy1 - energy1[0]) / abs(energy1[0]), label='Энергия Планеты 1', alpha=0.7)
axs[1].plot(t_3body, (energy2 - energy2[0]) / abs(energy2[0]), label='Энергия Планеты 2', alpha=0.7)
axs[1].plot(t_3body, (total_energy - total_energy[0]) / abs(total_energy[0]), label='Полная энергия системы',
            linewidth=2, color='black')
axs[1].set_title('Сохранение энергии в системе')
axs[1].set_xlabel('Время (годы)')
axs[1].set_ylabel('Относительное изменение E')
axs[1].grid(True)
axs[1].legend(loc='upper right')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
