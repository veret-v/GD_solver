#!/usr/bin/env python3
"""
Визуализация результатов расчёта вихря Тейлора–Грина.
Создаёт анимацию (GIF) и график затухания кинетической энергии.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import pandas as pd
import glob
import re
import argparse
import os

# ----------------------------------------------------------------------
# Аналитическое решение для затухания кинетической энергии
# ----------------------------------------------------------------------
def analytical_kinetic_energy(t, nu, L=2*np.pi):
    """
    Для вихря Тейлора–Грина в квадрате [0, 2π]²:
    E(t) = E0 * exp(-4 nu t)
    где E0 = 0.5 * (среднее по объёму от u²+v²) при t=0.
    При начальном поле u = sin x cos y, v = -cos x sin y:
    <u²> = <v²> = 1/4, поэтому E0 = 1/4.
    """
    E0 = 0.25
    return E0 * np.exp(-4.0 * nu * t)

# ----------------------------------------------------------------------
# Чтение CSV
# ----------------------------------------------------------------------
def read_csv(filepath):
    df = pd.read_csv(filepath)
    if not {'x', 'y', 'u', 'v'}.issubset(df.columns):
        raise ValueError(f"В файле {filepath} нет нужных столбцов x,y,u,v")
    return df

def build_grid(df):
    xs = np.sort(df['x'].unique())
    ys = np.sort(df['y'].unique())
    u_tab = df.pivot(index='y', columns='x', values='u').reindex(index=ys, columns=xs).to_numpy()
    v_tab = df.pivot(index='y', columns='x', values='v').reindex(index=ys, columns=xs).to_numpy()
    X, Y = np.meshgrid(xs, ys)
    return X, Y, u_tab, v_tab

def compute_kinetic_energy(df):
    """Вычисление средней кинетической энергии по объёму (среднее от 0.5*(u²+v²))."""
    u = df['u'].values
    v = df['v'].values
    return 0.5 * np.mean(u**2 + v**2)

# ----------------------------------------------------------------------
# Построение поля скорости
# ----------------------------------------------------------------------
def plot_field(ax, X, Y, U, V, title, vmax=None, add_quiver=True):
    speed = np.sqrt(U**2 + V**2)
    if vmax is None:
        vmax = np.nanmax(speed)
    cf = ax.contourf(X, Y, speed, levels=40, cmap='plasma', norm=Normalize(vmin=0, vmax=vmax))
    if add_quiver:
        step_x = max(1, X.shape[1] // 20)
        step_y = max(1, X.shape[0] // 20)
        ax.quiver(X[::step_y, ::step_x], Y[::step_y, ::step_x],
                  U[::step_y, ::step_x], V[::step_y, ::step_x],
                  color='white', scale=15, width=0.002)
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    return cf

# ----------------------------------------------------------------------
# Анимация
# ----------------------------------------------------------------------
def create_animation(files_info, output_file, fps=10):
    if not files_info:
        print("Нет файлов для анимации")
        return
    frames = []
    global_vmax = 0.0
    for fpath, tval in files_info:
        df = read_csv(fpath)
        X, Y, U, V = build_grid(df)
        speed = np.sqrt(U**2 + V**2)
        global_vmax = max(global_vmax, np.nanmax(speed))
        frames.append((X, Y, U, V, tval))
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    def update(frame_data):
        ax.clear()
        X, Y, U, V, tval = frame_data
        plot_field(ax, X, Y, U, V, f"t = {tval:.4f}", vmax=global_vmax, add_quiver=True)
        return []
    ani = animation.FuncAnimation(fig, update, frames=frames,
                                  interval=1000//fps, blit=False, repeat=True)
    ani.save(output_file, writer='pillow', fps=fps)
    plt.close(fig)
    print(f"GIF сохранён в {output_file}")

# ----------------------------------------------------------------------
# График сходимости (энергия)
# ----------------------------------------------------------------------
def plot_energy_convergence(times, E_num, nu, output_prefix):
    """
    Строит график max|u| / max|u0| от времени.
    Параметр E_num здесь больше не нужен и оставлен только ради совместимости сигнатуры.
    """
    if len(times) == 0:
        print("Нет данных для графика")
        return
    
    times_theory = np.linspace(0, max(times), 300)
    ratio_theory = np.exp(-2 * (np.pi**2) * nu * times_theory)
    # --- Конец добавления ---

    
    # Строим теоретическую кривую (пунктирная линия)
    

    times = np.asarray(times, dtype=float)

    # Предполагается, что E_num на самом деле передаётся как список max|u| по кадрам
    umax = np.asarray(E_num, dtype=float)

    if umax.size == 0:
        print("Нет данных для графика")
        return

    u0 = umax[0]
    if u0 == 0:
        raise ValueError("Начальное значение max|u0| равно 0, нормировка невозможна")

    ratio = umax / u0

    plt.figure(figsize=(8, 5))
    plt.plot(times_theory, ratio_theory, 'b--', label="theory: $e^{-2\pi^2 \\nu t}$")
    plt.plot(times, ratio, 'ro-', label="numeric")
    plt.xlabel('Время t')
    plt.ylabel(r"max |u| / max |u_0|")
    plt.title('Затухание максимальной скорости')
    plt.grid(True, alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_umax.png", dpi=150)
    plt.close()
    print(f"График max|u| сохранён в {output_prefix}_umax.png")

# ----------------------------------------------------------------------
# Точка входа
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Визуализация вихря Тейлора–Грина')
    parser.add_argument('--dir', type=str, default='taylor_green',
                        help='Директория с CSV-файлами')
    parser.add_argument('--pattern', type=str, default='step_*.csv',
                        help='Шаблон имён файлов')
    parser.add_argument('--output', type=str, default='taylor_green_animation.gif',
                        help='Имя выходного GIF')
    parser.add_argument('--fps', type=int, default=10,
                        help='Кадров в секунду')
    parser.add_argument('--nu', type=float, default=0.025,
                        help='Кинематическая вязкость (для аналитики)')
    args = parser.parse_args()

    pattern = os.path.join(args.dir, args.pattern)
    files = glob.glob(pattern)
    if not files:
        print(f"Файлы не найдены: {pattern}")
        return

    def extract_step(fname):
        m = re.search(r'step_(\d+)', fname)
        return int(m.group(1)) if m else 0
    files.sort(key=extract_step)

    file_info = []
    times = []
    umax_values = []

    for f in files:
        df = read_csv(f)
        t_match = re.search(r'_t_([\d\.]+)', f)
        t_val = float(t_match.group(1)[:-1]) if t_match else 0.0
        file_info.append((f, t_val))
        times.append(t_val)
        umax_values.append(np.max(np.abs(df['u'].values)))

    if not file_info:
        print("Нет данных для визуализации")
        return

    create_animation(file_info, args.output, args.fps)
    plot_energy_convergence(
        np.array(times),
        np.array(umax_values),
        args.nu,
        os.path.splitext(args.output)[0]
    )

if __name__ == '__main__':
    main()