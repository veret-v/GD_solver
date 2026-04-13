#!/usr/bin/env python3
"""
Визуализация результатов расчёта каверны (SIMPLE/PISO/PIMPLE).
Создаёт анимацию (GIF) и графики сравнения с Ghia et al. (1982).
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
# Данные Ghia et al. (Re=100)
# ----------------------------------------------------------------------
GHIA_Y_U = np.array([0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
                     0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
                     0.9688, 0.9766, 1.0000])
GHIA_U   = np.array([0.00000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150,
                     -0.15662, -0.21090, -0.20581, -0.13641,  0.00332,  0.23151,
                      0.68717,  0.73722,  0.78871,  0.84123,  1.00000])

GHIA_X_V = np.array([0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
                     0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
                     0.9609, 0.9688, 1.0000])
GHIA_V   = np.array([0.00000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077,
                     0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914,
                     -0.10313, -0.08864, -0.07391, -0.05906, 0.00000])

# ----------------------------------------------------------------------
# Чтение CSV-файла
# ----------------------------------------------------------------------
def read_csv_1(filepath):
    data = pd.read_csv(filepath)
    if data.ndim == 1:
        return None, None, None, None, None
    x = np.array(data['x'])
    y = np.array(data['y'])
    u = np.array(data['u'])
    v = np.array(data['v'])
    p = np.array(data['p'])
    return x, y, u, v, p

def read_csv(filepath):
    df = pd.read_csv(filepath)
    if not {'x', 'y', 'u', 'v'}.issubset(df.columns):
        raise ValueError(f"В файле {filepath} нет нужных столбцов x,y,u,v")
    return df

def build_grid(df):
    xs = np.sort(df['x'].unique())
    ys = np.sort(df['y'].unique())

    u_tab = (
        df.pivot(index='y', columns='x', values='u')
        .reindex(index=ys, columns=xs)
        .to_numpy()
    )
    v_tab = (
        df.pivot(index='y', columns='x', values='v')
        .reindex(index=ys, columns=xs)
        .to_numpy()
    )

    X, Y = np.meshgrid(xs, ys)
    return X, Y, u_tab, v_tab

 
def plot_field(ax, X, Y, U, V, title, vmax=None, add_quiver=True):
    speed = np.sqrt(U**2 + V**2)

    if vmax is None:
        vmax = np.nanmax(speed)

    cf = ax.contourf(
        X, Y, speed,
        levels=40,
        cmap='plasma',
        norm=Normalize(vmin=0, vmax=vmax)
    )

    if add_quiver:
        step_x = max(1, X.shape[1] // 16)
        step_y = max(1, X.shape[0] // 16)

        ax.quiver(
            X[::step_y, ::step_x],
            Y[::step_y, ::step_x],
            U[::step_y, ::step_x],
            V[::step_y, ::step_x],
            color='white',
            scale=15,
            width=0.002
        )

    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    return cf


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

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=1000 // fps,
        blit=False,
        repeat=True
    )

    ani.save(output_file, writer='pillow', fps=fps)
    plt.close(fig)
    print(f"GIF сохранён в {output_file}")

# ----------------------------------------------------------------------
# Сравнение с Ghia et al.
# ----------------------------------------------------------------------
def plot_ghia_comparison(X, Y, U, V, output_prefix):
    x_center = 0.5
    y_center = 0.5

    # Индексы ближайших точек
    print(np.abs(X - x_center))
    ix = X[np.argmin(np.abs(X - x_center))]
    jy = Y[np.argmin(np.abs(Y - y_center))]
    
    print(ix, jy, sep=" ")

    u_profile = U[np.abs(X - ix) < 1e-2]          # u вдоль вертикальной линии x=0.5
    y_coord   = Y[np.abs(X - ix) < 1e-2]

    v_profile = V[np.abs(Y - jy) < 1e-2]          # v вдоль горизонтальной линии y=0.5
    x_coord   = X[np.abs(Y - jy) < 1e-2]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # u(y)
    axes[0].plot(u_profile, y_coord, 'r-', label='Расчёт')
    axes[0].plot(GHIA_U, GHIA_Y_U, 'ko', markersize=4, label='Ghia et al.')
    axes[0].set_xlabel('u')
    axes[0].set_ylabel('y')
    axes[0].set_title('u(y) на x=0.5')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # v(x)
    axes[1].plot(x_coord, v_profile, 'r-', label='Расчёт')
    axes[1].plot(GHIA_X_V, GHIA_V, 'ko', markersize=4, label='Ghia et al.')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('v')
    axes[1].set_title('v(x) на y=0.5')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(f"{output_prefix}_ghia_comparison.png", dpi=150)
    plt.close(fig)
    print(f"Графики сравнения сохранены в {output_prefix}_ghia_comparison.png")

# ----------------------------------------------------------------------
# Точка входа
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Визуализация каверны')
    parser.add_argument('--dir', type=str, default='steps',
                        help='Директория с CSV-файлами')
    parser.add_argument('--pattern', type=str, default='step_*.csv',
                        help='Шаблон имён файлов')
    parser.add_argument('--output', type=str, default='cavity_animation.gif',
                        help='Имя выходного GIF')
    parser.add_argument('--fps', type=int, default=10,
                        help='Кадров в секунду')
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
    last_data = None
    for f in files:
        x, y, u, v, _ = read_csv_1(f)
        if x is None:
            continue
        t_match = re.search(r'_t_([\d\.]+)', f)
        t_val = float(t_match.group(1)[:-1]) if t_match else 0.0
        file_info.append((f, t_val))
        last_data = (x, y, u, v)

    if not file_info:
        print("Нет данных для анимации")
        return

    create_animation(file_info, args.output, args.fps)

    if last_data:
        x, y, u, v = last_data
        prefix = os.path.splitext(args.output)[0]
        plot_ghia_comparison(x, y, u, v, prefix)

if __name__ == '__main__':
    main()