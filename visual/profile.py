import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import imageio
import os
import re

def create_animation_2D(filenames, output_gif='animation_2d.gif', gamma=1.4):
    """
    Создаёт GIF-анимацию из серии 2D CSV-файлов с результатами расчёта.
    
    Параметры
    ----------
    filenames : list
        Список путей к CSV-файлам (каждый файл – один временной слой).
    output_gif : str
        Имя выходного GIF-файла.
    gamma : float
        Показатель адиабаты для вычисления внутренней энергии.
    """
    # Извлекаем время из имён файлов и сортируем
    data_by_time = {}
    for filename in filenames:
        match = re.search(r'time_([\d.]+)\.csv', filename)
        if match:
            time = float(match.group(1))
            df = pd.read_csv(filename)
            data_by_time[time] = df
        else:
            print(f"Предупреждение: не удалось извлечь время из имени {filename}, пропускаем.")

    if not data_by_time:
        print("Не найдено ни одного файла с корректной временной меткой.")
        return

    times = sorted(data_by_time.keys())
    print(f"Найдено временных слоёв: {len(times)}")

    # Определим уникальные координаты сетки по первому файлу (предполагаем, что сетка постоянна)
    sample_df = data_by_time[times[0]]
    x_vals = np.sort(sample_df['x'].unique())
    y_vals = np.sort(sample_df['y'].unique())
    nx, ny = len(x_vals), len(y_vals)
    print(f"Сетка: {nx} x {ny} (x от {x_vals.min():.3f} до {x_vals.max():.3f}, "
          f"y от {y_vals.min():.3f} до {y_vals.max():.3f})")

    # Подготовка к сохранению кадров
    frame_filenames = []
    for i, t in enumerate(times):
        df = data_by_time[t]

        # Преобразуем табличные данные в 2D-массивы для каждой величины
        try:
            rho_num = df.pivot(index='y', columns='x', values='rho').values
            rho_ex = df.pivot(index='y', columns='x', values='rho_exact').values
            p_num = df.pivot(index='y', columns='x', values='p').values
            p_ex = df.pivot(index='y', columns='x', values='p_exact').values
            u_num = df.pivot(index='y', columns='x', values='u').values
            u_ex = df.pivot(index='y', columns='x', values='u_exact').values
            v_num = df.pivot(index='y', columns='x', values='v').values
            v_ex = df.pivot(index='y', columns='x', values='v_exact').values
        except KeyError as e:
            print(f"Ошибка: в файле для времени {t} отсутствует столбец {e}")
            continue

        # Вычисляем внутреннюю энергию
        e_num = p_num / ((gamma - 1) * rho_num)
        e_ex = p_ex / ((gamma - 1) * rho_ex)

        # Создаём сетку для отрисовки (матрицы x и y)
        X, Y = np.meshgrid(x_vals, y_vals)

        # Строим фигуру с 4 подграфиками (можно настроить под свои нужды)
        fig, axs = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'2D поля в момент времени t = {t:.6f} с', fontsize=16)

        # Настройка пределов для цветовых шкал (можно брать глобальные min/max по всем временам)
        # Для простоты используем локальные min/max каждого кадра
        plots = [
            (rho_num, rho_ex, r'Плотность $\rho$', axs[0, 0]),
            (p_num,   p_ex,   r'Давление $p$',      axs[0, 1]),
            (u_num,   u_ex,   r'Скорость $u$',      axs[1, 0]),
            (v_num,   v_ex,   r'Скорость $v$',      axs[1, 1])
        ]

        for num, ex, title, ax in plots:
            # Рисуем численное решение (заливка)
            im = ax.pcolormesh(X, Y, num, shading='auto', cmap='viridis')
            # Можно добавить контуры точного решения (опционально)
            # ax.contour(X, Y, ex, colors='red', linewidths=0.5, linestyles='--')
            ax.set_title(title)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_aspect('equal')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # Альтернативно: можно показывать ошибку или другие комбинации
        # Здесь оставлено для простоты – численные поля.

        plt.tight_layout()

        # Сохраняем кадр
        frame_name = f'frame_{i:03d}.png'
        plt.savefig(frame_name, dpi=100, bbox_inches='tight')
        plt.close()
        frame_filenames.append(frame_name)
        print(f"Кадр {i+1}/{len(times)} сохранён (t = {t:.6f})")

    # Создаём GIF
    if frame_filenames:
        print("Создание GIF...")
        with imageio.get_writer(output_gif, mode='I', duration=0.2) as writer:
            for fname in frame_filenames:
                image = imageio.imread(fname)
                writer.append_data(image)
        # Удаляем временные файлы
        for fname in frame_filenames:
            os.remove(fname)
        print(f"GIF успешно создан: {output_gif}")
    else:
        print("Не создано ни одного кадра.")


if __name__ == "__main__":
    # Путь к папке с файлами (измените на свой)
    path = r'C:\solvver\GD_solver\output\steps'

    if not os.path.exists(path):
        print(f"Ошибка: путь {path} не существует!")
        exit()

    # Ищем все CSV-файлы
    filenames = glob.glob(os.path.join(path, "*.csv"))
    if not filenames:
        # Пробуем текущую папку
        filenames = glob.glob("step_*_time_*.csv")
        if not filenames:
            print("CSV-файлы не найдены.")
            exit()

    print(f"Найдено файлов: {len(filenames)}")
    # Для проверки покажем структуру первого
    sample_df = pd.read_csv(filenames[0])
    print("Столбцы первого файла:", sample_df.columns.tolist())
    print("Количество строк:", len(sample_df))

    create_animation_2D(filenames, output_gif='animation_2d.gif', gamma=1.4)