import numpy as np
import matplotlib.pyplot as plt

# Параметры системы
q = 1.0      # Заряд
B = 1.0      # Магнитное поле
m = 1.0      # Масса
omega = q * B / m  # Циклотронная частота

# Начальные условия
x0, y0 = 0.0, 0.0
vx0, vy0 = 0.0, 10.0

# Время и шаг интегрирования
t_start, t_end = 0.0, 10.0
dt = 0.01
steps = int((t_end - t_start) / dt)

# Массивы для результатов
t_values = np.linspace(t_start, t_end, steps)
x_values = np.zeros(steps)
y_values = np.zeros(steps)
vx_values = np.zeros(steps)
vy_values = np.zeros(steps)

# Начальные значения
x_values[0], y_values[0] = x0, y0
vx_values[0], vy_values[0] = vx0, vy0

# Функция правых частей системы
def f(t, u):
    x, y, vx, vy = u
    dxdt = vx
    dydt = vy
    dvxdt = omega * vy
    dvydt = -omega * vx
    return np.array([dxdt, dydt, dvxdt, dvydt])

# Интегрирование методом RK4
for i in range(steps - 1):
    t = t_values[i]
    u_current = np.array([x_values[i], y_values[i], vx_values[i], vy_values[i]])

    k1 = f(t, u_current) * dt
    k2 = f(t + dt/2, u_current + k1/2) * dt
    k3 = f(t + dt/2, u_current + k2/2) * dt
    k4 = f(t + dt, u_current + k3) * dt

    u_next = u_current + (k1 + 2*k2 + 2*k3 + k4) / 6
    x_values[i+1], y_values[i+1], vx_values[i+1], vy_values[i+1] = u_next

# Визуализация траектории
plt.figure(figsize=(8, 6))
plt.plot(x_values, y_values, label='Траектория частицы')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Движение заряженной частицы в магнитном поле')
plt.grid(True)
plt.axis('equal')
plt.legend()
plt.show()
