import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import glob
import imageio.v2 as imageio
import os
import re
import traceback

def load_method_tag(project_root):
    equation_type = -1
    cfg_path = os.path.join(project_root, "configs", "input.ini")
    if os.path.exists(cfg_path):
        in_system = False
        with open(cfg_path, "r", encoding="utf-8", errors="ignore") as f:
            for raw in f:
                line = raw.split("#", 1)[0].strip()
                if not line:
                    continue
                if line.startswith("[") and line.endswith("]"):
                    in_system = (line[1:-1].strip().lower() == "system")
                    continue
                if not in_system or "=" not in line:
                    continue
                key, value = [s.strip() for s in line.split("=", 1)]
                if key.lower() == "equation_type":
                    try:
                        equation_type = int(value)
                    except ValueError:
                        pass
    method_names = {
        0: "Godunov", 1: "Kolgan", 2: "Rodionov", 3: "HLL",
        4: "HLLC", 5: "Rusanov", 6: "Osher", 7: "Roe", 8: "FLIC", 9: "Mader2DE",
    }
    return f"eq={equation_type} ({method_names.get(equation_type, 'Unknown')})"


def has_valid_exact(df):
    exact_cols = ["rho_exact", "u_exact", "v_exact", "p_exact"]
    if not all(col in df.columns for col in exact_cols):
        return False
    return np.isfinite(df[exact_cols].to_numpy(dtype=float)).any()


def build_mask(df, x_vals, y_vals):
    if 'is_solid' not in df.columns:
        return np.zeros((len(y_vals), len(x_vals)), dtype=bool)
    pivot = df.pivot(index='y', columns='x', values='is_solid')
    pivot = pivot.reindex(index=y_vals, columns=x_vals)
    return pivot.values > 0.5


def get_2d(df, x_vals, y_vals, var, solid_mask):
    pivot = df.pivot(index='y', columns='x', values=var)
    pivot = pivot.reindex(index=y_vals, columns=x_vals)
    arr = pivot.values.astype(float)
    mask = solid_mask | ~np.isfinite(arr)
    return np.ma.array(arr, mask=mask)


def create_2d_animation(filenames, method_tag, output_gif='heatmap_2d.gif'):
    """
    Создаёт GIF-анимацию, показывающую 2D-распределения величин.
    """
    data = {}
    for filename in filenames:
        match = re.search(r'time_([\d.]+)\.csv', filename)
        if match:
            try:
                time = float(match.group(1))
                if not np.isfinite(time):
                    print(f"Предупреждение: время {time} не конечное, пропускаем {filename}")
                    continue
                df = pd.read_csv(filename)
                data[time] = df
            except Exception as e:
                print(f"Ошибка при чтении {filename}: {e}")
                continue

    times = sorted(data.keys())
    if not times:
        print("Нет корректных данных по времени!")
        return

    frame_files = []
    for i, t in enumerate(times):
        try:
            df = data[t]
            x_vals = np.sort(df['x'].unique())
            y_vals = np.sort(df['y'].unique())
            X, Y = np.meshgrid(x_vals, y_vals)
            solid_mask = build_mask(df, x_vals, y_vals)

            # Загружаем данные
            rho = get_2d(df, x_vals, y_vals, 'rho', solid_mask)
            u = get_2d(df, x_vals, y_vals, 'u', solid_mask)
            v = get_2d(df, x_vals, y_vals, 'v', solid_mask)
            p = get_2d(df, x_vals, y_vals, 'p', solid_mask)
            exact_available = has_valid_exact(df)
            if exact_available:
                rho_ex = get_2d(df, x_vals, y_vals, 'rho_exact', solid_mask)
                u_ex = get_2d(df, x_vals, y_vals, 'u_exact', solid_mask)
                v_ex = get_2d(df, x_vals, y_vals, 'v_exact', solid_mask)
                p_ex = get_2d(df, x_vals, y_vals, 'p_exact', solid_mask)

            gamma = 1.4
            e = p / (rho * (gamma - 1))
            if exact_available:
                e_ex = p_ex / (rho_ex * (gamma - 1))

            # Создаём фигуру
            if exact_available:
                fig, axs = plt.subplots(3, 2, figsize=(14, 18))
            else:
                fig, axs = plt.subplots(3, 1, figsize=(10, 18))
            # Упрощённый заголовок (без форматирования, чтобы избежать возможной проблемы)
            fig.suptitle(f'Time: {t} s | {method_tag}', fontsize=16)

            if exact_available:
                # Плотность
                vmin_rho = float(min(np.ma.min(rho), np.ma.min(rho_ex)))
                vmax_rho = float(max(np.ma.max(rho), np.ma.max(rho_ex)))
                im1 = axs[0,0].pcolormesh(X, Y, rho, shading='auto', norm=Normalize(vmin_rho, vmax_rho), cmap='viridis')
                axs[0,0].set_title('Density (numerical)')
                axs[0,0].set_xlabel('x'); axs[0,0].set_ylabel('y')
                plt.colorbar(im1, ax=axs[0,0])

                im2 = axs[0,1].pcolormesh(X, Y, rho_ex, shading='auto', norm=Normalize(vmin_rho, vmax_rho), cmap='viridis')
                axs[0,1].set_title('Density (exact)')
                axs[0,1].set_xlabel('x'); axs[0,1].set_ylabel('y')
                plt.colorbar(im2, ax=axs[0,1])

                # x-скорость
                vmin_u = float(min(np.ma.min(u), np.ma.min(u_ex)))
                vmax_u = float(max(np.ma.max(u), np.ma.max(u_ex)))
                im3 = axs[1,0].pcolormesh(X, Y, u, shading='auto', norm=Normalize(vmin_u, vmax_u), cmap='plasma')
                axs[1,0].set_title('u-velocity (numerical)')
                axs[1,0].set_xlabel('x'); axs[1,0].set_ylabel('y')
                plt.colorbar(im3, ax=axs[1,0])

                im4 = axs[1,1].pcolormesh(X, Y, u_ex, shading='auto', norm=Normalize(vmin_u, vmax_u), cmap='plasma')
                axs[1,1].set_title('u-velocity (exact)')
                axs[1,1].set_xlabel('x'); axs[1,1].set_ylabel('y')
                plt.colorbar(im4, ax=axs[1,1])

                # Внутренняя энергия
                vmin_e = float(min(np.ma.min(e), np.ma.min(e_ex)))
                vmax_e = float(max(np.ma.max(e), np.ma.max(e_ex)))
                im5 = axs[2,0].pcolormesh(X, Y, e, shading='auto', norm=Normalize(vmin_e, vmax_e), cmap='inferno')
                axs[2,0].set_title('Internal energy (numerical)')
                axs[2,0].set_xlabel('x'); axs[2,0].set_ylabel('y')
                plt.colorbar(im5, ax=axs[2,0])

                im6 = axs[2,1].pcolormesh(X, Y, e_ex, shading='auto', norm=Normalize(vmin_e, vmax_e), cmap='inferno')
                axs[2,1].set_title('Internal energy (exact)')
                axs[2,1].set_xlabel('x'); axs[2,1].set_ylabel('y')
                plt.colorbar(im6, ax=axs[2,1])
            else:
                fields = [
                    (rho, 'Density (numerical)', 'viridis'),
                    (u, 'u-velocity (numerical)', 'plasma'),
                    (e, 'Internal energy (numerical)', 'inferno'),
                ]
                for ax, (arr, title, cmap) in zip(axs, fields):
                    im = ax.pcolormesh(X, Y, arr, shading='auto', cmap=cmap)
                    ax.set_title(title)
                    ax.set_xlabel('x')
                    ax.set_ylabel('y')
                    plt.colorbar(im, ax=ax)

            plt.tight_layout()
            frame_name = f'frame_2d_{i:03d}.png'
            plt.savefig(frame_name, dpi=100, bbox_inches='tight')
            plt.close()
            frame_files.append(frame_name)
            print(f"Кадр {i+1}/{len(times)} (время {t}) успешно создан")

        except Exception as e:
            print(f"Ошибка при обработке кадра {i} (время {t}):")
            traceback.print_exc()
            # Продолжаем со следующим кадром, не останавливая цикл

    if frame_files:
        print("Создание GIF...")
        try:
            with imageio.get_writer(output_gif, mode='I', duration=0.5) as writer:
                for fname in frame_files:
                    image = imageio.imread(fname)
                    writer.append_data(image)
            for fname in frame_files:
                os.remove(fname)
            print(f"GIF сохранён: {output_gif}")
        except Exception as e:
            print("Ошибка при создании GIF:")
            traceback.print_exc()
    else:
        print("Не создано ни одного кадра!")

if __name__ == "__main__":
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    method_tag = load_method_tag(project_root)
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "output", "steps"))   # измените при необходимости
    if not os.path.exists(path):
        print(f"Путь {path} не найден, ищем файлы в текущей директории...")
        filenames = glob.glob("step_*_time_*.csv")
    else:
        filenames = glob.glob(os.path.join(path, "step_*_time_*.csv"))

    if not filenames:
        print("CSV-файлы не найдены!")
        exit()

    print(f"Найдено {len(filenames)} файлов")
    create_2d_animation(filenames, method_tag, output_gif='heatmap_2d.gif')
