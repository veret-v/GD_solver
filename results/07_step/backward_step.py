#!/usr/bin/env python3
"""
Визуализация результатов расчёта обратной ступеньки.
Создаёт анимацию (GIF) — тепловая карта скорости в выбранной области
с вырезом на месте ступеньки, и профили u(y) в нескольких сечениях x = const.
"""

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
from matplotlib.colors import Normalize, PowerNorm


def read_csv(filepath):
    df = pd.read_csv(filepath)
    required = {'x', 'y', 'u', 'v'}
    if not required.issubset(df.columns):
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


def parse_time_from_name(fname):
    m = re.search(r'_t_([\d\.]+)', fname)
    if not m:
        return 0.0
    s = m.group(1)
    if s.endswith('.'):
        s = s[:-1]
    return float(s)


def extract_step_number(fname):
    m = re.search(r'step_(\d+)', os.path.basename(fname))
    return int(m.group(1)) if m else 0


def make_step_mask(X, Y, step_x, step_y):
    return (X <= step_x) & (Y <= step_y)


def apply_step_mask(Z, mask):
    Zm = np.ma.array(Z, copy=True)
    Zm[mask] = np.ma.masked
    return Zm


def add_step_patch(ax, step_x, step_y):
    rect = plt.Rectangle(
        (0.0, 0.0),
        step_x,
        step_y,
        facecolor='white',
        edgecolor='black',
        linewidth=1.0,
        zorder=10
    )
    ax.add_patch(rect)


def plot_heatmap(ax, X, Y, U, V, title, vmax=None, xlim=None, ylim=None,
                 step_x=1.0, step_y=0.5, show_step=True):
    speed = np.sqrt(U**2 + V**2)

    if vmax is None:
        vmax = np.nanpercentile(speed, 98)

    mask = make_step_mask(X, Y, step_x, step_y)
    speed_masked = apply_step_mask(speed, mask)

    # Более яркая колормэпа с уменьшенным gamma и яркими пределами
    cf = ax.contourf(
        X, Y, speed_masked,
        levels=80,
        cmap='plasma',  # Поменял на plasma - более яркая и контрастная
        norm=PowerNorm(gamma=0.4, vmin=0, vmax=vmax),  # Уменьшил gamma до 0.4
        corner_mask=False,
        extend='max'  # Добавляет яркий цвет для значений выше vmax
    )

    if show_step:
        add_step_patch(ax, step_x, step_y)

    if xlim is not None:
        ax.set_xlim(*xlim)
    else:
        ax.set_xlim(X.min(), X.max())

    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        ax.set_ylim(Y.min(), Y.max())

    ax.set_aspect('equal', adjustable='box')
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    
    # Добавляем colorba
    
    return cf


def create_animation(files_info, output_file, fps=10, xlim=None, ylim=None,
                     step_x=1.0, step_y=0.5):
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

    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)

    def update(frame_data):
        ax.clear()
        X, Y, U, V, tval = frame_data
        plot_heatmap(
            ax, X, Y, U, V,
            title=f"t = {tval:.4f}",
            vmax=global_vmax,
            xlim=xlim,
            ylim=ylim,
            step_x=step_x,
            step_y=step_y,
            show_step=True
        )
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


def plot_u_profiles(df, x_positions, output_prefix):
    xs = np.sort(df['x'].unique())
    ys = np.sort(df['y'].unique())

    u_tab = (
        df.pivot(index='y', columns='x', values='u')
        .reindex(index=ys, columns=xs)
        .to_numpy()
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, len(x_positions)))

    for x_target, color in zip(x_positions, colors):
        idx = int(np.argmin(np.abs(xs - x_target)))
        x_actual = xs[idx]
        u_profile = u_tab[:, idx]
        ax.plot(u_profile, ys, color=color, label=f'x = {x_actual:.3f}')

    ax.set_xlabel('u')
    ax.set_ylabel('y')
    ax.set_title('Профили горизонтальной скорости u(y)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{output_prefix}_u_profiles.png", dpi=150)
    plt.close(fig)
    print(f"Профили скорости сохранены в {output_prefix}_u_profiles.png")


def main():
    parser = argparse.ArgumentParser(description='Визуализация обратной ступеньки')
    parser.add_argument('--dir', type=str, default='backward_step',
                        help='Директория с CSV-файлами')
    parser.add_argument('--pattern', type=str, default='step_*.csv',
                        help='Шаблон имён файлов')
    parser.add_argument('--output', type=str, default='backward_step_animation.gif',
                        help='Имя выходного GIF')
    parser.add_argument('--fps', type=int, default=10,
                        help='Кадров в секунду')
    parser.add_argument('--x_sections', type=float, nargs='+',
                        default=[4.0, 5.0, 7.0, 10.0, 14.0],
                        help='Координаты x для построения профилей u(y)')
    parser.add_argument('--xlim', type=float, nargs=2, default=None,
                        help='Границы по x для выреза на GIF (xmin xmax)')
    parser.add_argument('--ylim', type=float, nargs=2, default=None,
                        help='Границы по y для выреза на GIF (ymin ymax)')
    parser.add_argument('--step_x', type=float, default=2.5,
                        help='Длина ступеньки по x')
    parser.add_argument('--step_y', type=float, default=0.5,
                        help='Высота ступеньки по y')
    args = parser.parse_args()

    pattern = os.path.join(args.dir, args.pattern)
    files = glob.glob(pattern)
    if not files:
        print(f"Файлы не найдены: {pattern}")
        return

    files.sort(key=extract_step_number)

    file_info = []
    last_df = None

    for f in files:
        df = read_csv(f)
        t_val = parse_time_from_name(f)
        file_info.append((f, t_val))
        last_df = df

    if not file_info:
        print("Нет данных для визуализации")
        return

    xlim = tuple(args.xlim) if args.xlim else None
    ylim = tuple(args.ylim) if args.ylim else None

    create_animation(
        file_info,
        args.output,
        args.fps,
        xlim=xlim,
        ylim=ylim,
        step_x=args.step_x,
        step_y=args.step_y
    )

    if last_df is not None:
        plot_u_profiles(last_df, args.x_sections, os.path.splitext(args.output)[0])


if __name__ == '__main__':
    main()