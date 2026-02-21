import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import glob
import imageio
import os
import re
import traceback

def create_2d_animation(filenames, output_gif='animation_2d.gif'):
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

            def get_2d(var):
                pivot = df.pivot(index='y', columns='x', values=var)
                pivot = pivot.reindex(index=y_vals, columns=x_vals)
                return pivot.values

            # Загружаем данные
            rho = get_2d('rho')
            rho_ex = get_2d('rho_exact')
            u = get_2d('u')
            u_ex = get_2d('u_exact')
            v = get_2d('v')
            v_ex = get_2d('v_exact')
            p = get_2d('p')
            p_ex = get_2d('p_exact')

            # Проверка на наличие нечисловых значений
            for name, arr in [('rho', rho), ('rho_ex', rho_ex), ('u', u), ('u_ex', u_ex),
                              ('v', v), ('v_ex', v_ex), ('p', p), ('p_ex', p_ex)]:
                if not np.all(np.isfinite(arr)):
                    print(f"Кадр {i} (время {t}): массив {name} содержит нечисловые значения, пропускаем")
                    raise ValueError("Non-finite data")

            gamma = 1.4
            e = p / (rho * (gamma - 1))
            e_ex = p_ex / (rho_ex * (gamma - 1))
            if not (np.all(np.isfinite(e)) and np.all(np.isfinite(e_ex))):
                print(f"Кадр {i}: внутренняя энергия содержит inf/nan, пропускаем")
                raise ValueError("Non-finite internal energy")

            # Создаём фигуру
            fig, axs = plt.subplots(3, 2, figsize=(14, 18))
            # Упрощённый заголовок (без форматирования, чтобы избежать возможной проблемы)
            fig.suptitle(f'Time: {t} s', fontsize=16)

            # Плотность
            vmin_rho = min(rho.min(), rho_ex.min())
            vmax_rho = max(rho.max(), rho_ex.max())
            im1 = axs[0,0].pcolormesh(X, Y, rho, shading='auto', norm=Normalize(vmin_rho, vmax_rho), cmap='viridis')
            axs[0,0].set_title('Density (numerical)')
            axs[0,0].set_xlabel('x'); axs[0,0].set_ylabel('y')
            plt.colorbar(im1, ax=axs[0,0])

            im2 = axs[0,1].pcolormesh(X, Y, rho_ex, shading='auto', norm=Normalize(vmin_rho, vmax_rho), cmap='viridis')
            axs[0,1].set_title('Density (exact)')
            axs[0,1].set_xlabel('x'); axs[0,1].set_ylabel('y')
            plt.colorbar(im2, ax=axs[0,1])

            # x-скорость
            vmin_u = min(u.min(), u_ex.min())
            vmax_u = max(u.max(), u_ex.max())
            im3 = axs[1,0].pcolormesh(X, Y, u, shading='auto', norm=Normalize(vmin_u, vmax_u), cmap='plasma')
            axs[1,0].set_title('u-velocity (numerical)')
            axs[1,0].set_xlabel('x'); axs[1,0].set_ylabel('y')
            plt.colorbar(im3, ax=axs[1,0])

            im4 = axs[1,1].pcolormesh(X, Y, u_ex, shading='auto', norm=Normalize(vmin_u, vmax_u), cmap='plasma')
            axs[1,1].set_title('u-velocity (exact)')
            axs[1,1].set_xlabel('x'); axs[1,1].set_ylabel('y')
            plt.colorbar(im4, ax=axs[1,1])

            # Внутренняя энергия
            vmin_e = min(e.min(), e_ex.min())
            vmax_e = max(e.max(), e_ex.max())
            im5 = axs[2,0].pcolormesh(X, Y, e, shading='auto', norm=Normalize(vmin_e, vmax_e), cmap='inferno')
            axs[2,0].set_title('Internal energy (numerical)')
            axs[2,0].set_xlabel('x'); axs[2,0].set_ylabel('y')
            plt.colorbar(im5, ax=axs[2,0])

            im6 = axs[2,1].pcolormesh(X, Y, e_ex, shading='auto', norm=Normalize(vmin_e, vmax_e), cmap='inferno')
            axs[2,1].set_title('Internal energy (exact)')
            axs[2,1].set_xlabel('x'); axs[2,1].set_ylabel('y')
            plt.colorbar(im6, ax=axs[2,1])

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
    path = r'C:\solvver\GD_solver\output\steps'   # измените при необходимости
    if not os.path.exists(path):
        print(f"Путь {path} не найден, ищем файлы в текущей директории...")
        filenames = glob.glob("step_*_time_*.csv")
    else:
        filenames = glob.glob(os.path.join(path, "step_*_time_*.csv"))

    if not filenames:
        print("CSV-файлы не найдены!")
        exit()

    print(f"Найдено {len(filenames)} файлов")
    create_2d_animation(filenames, output_gif='animation_2d.gif')