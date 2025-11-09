import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import imageio
import os

path = '/home/jupyter-eds010'  # Вставить путь со всеми расчетами
filenames = glob.glob(path + "/*.csv")

# Находит время из названия файла. Пример: "solution_0.1.csv". Время 0.1
time = sorted([float(f.split('_')[-1].strip('.csv')) for f in filenames])
# Закидываем все файлы в один массив
dfs = [pd.read_csv(filename, sep=';') for filename in filenames]
# Словарь, где ключ - это время, а значения - это массив с параметрами (может быть неудобно)
data = {t: param for t, param in zip(time, dfs)}


def create_animation_1D(filenames, output_gif='animation.gif'):

    # Сортируем временные точки
    times = sorted(data.keys())

    filenames_write = []

    # Создаем кадры для каждого момента времени
    for i, t in enumerate(times):
        df = data[t]

        # Предполагаем структуру DataFrame: [x, velocity, density, pressure]
        x_coord = df.iloc[:, 0]  # первая колонка - координата
        velocity = df.iloc[:, 1]  # вторая колонка - скорость
        density = df.iloc[:, 2]  # третья колонка - плотность
        pressure = df.iloc[:, 3]  # четвертая колонка - давление

        fig, axs = plt.subplots(3, 1, figsize=(8, 10))
        fig.suptitle(f'Time: {t:.2f} s')

        # График скорости от координаты
        axs[0].plot(x_coord, velocity, 'b-', label='Velocity', linewidth=2)
        axs[0].set_ylabel('Velocity')
        axs[0].grid(True, alpha=0.3)
        axs[0].legend()

        # График плотности от координаты
        axs[1].plot(x_coord, density, 'g-', label='Density', linewidth=2)
        axs[1].set_ylabel('Density')
        axs[1].grid(True, alpha=0.3)
        axs[1].legend()

        # График давления от координаты
        axs[2].plot(x_coord, pressure, 'r-', label='Pressure', linewidth=2)
        axs[2].set_ylabel('Pressure')
        axs[2].set_xlabel('Coordinate')
        axs[2].grid(True, alpha=0.3)
        axs[2].legend()

        plt.tight_layout()

        # Сохраняем кадр
        filename = f'frame_{i:03d}.png'
        plt.savefig(filename, dpi=100, bbox_inches='tight')
        plt.close()
        filenames_write.append(filename)

    # Создаем GIF из кадров
    with imageio.get_writer(output_gif, mode='I', duration=0.5) as writer:
        for filename in filenames_write:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Удаляем временные файлы кадров
    for filename in filenames_write:
        os.remove(filename)


def create_animation_2D(filenames, output_gif='animation_2D.gif'):
    # Извлекаем время из названий файлов и создаем словарь с данными
    times = sorted([float(f.split('_')[-1].strip('.csv')) for f in filenames])
    dfs = [pd.read_csv(filename, sep=';') for filename in filenames]
    data = {t: df for t, df in zip(times, dfs)}

    # Сортируем временные точки
    times = sorted(data.keys())

    filenames_write = []

    # Создаем кадры для каждого момента времени
    for i, t in enumerate(times):
        df = data[t]

        # Предполагаем структуру DataFrame: [x, y, vx, vy, density, pressure]
        x_coord = df.iloc[:, 0]  # первая колонка - x координата
        y_coord = df.iloc[:, 1]  # вторая колонка - y координата
        vx = df.iloc[:, 2]       # третья колонка - скорость по x
        vy = df.iloc[:, 3]       # четвертая колонка - скорость по y
        density = df.iloc[:, 4]  # пятая колонка - плотность
        pressure = df.iloc[:, 5]  # шестая колонка - давление

        # Создаем сетку для тепловых карт
        # Находим уникальные значения координат
        x_unique = np.sort(x_coord.unique())
        y_unique = np.sort(y_coord.unique())

        # Создаем сетку координат
        X, Y = np.meshgrid(x_unique, y_unique)

        # Преобразуем данные в 2D массивы для тепловых карт
        # Для этого создаем pivot таблицы
        df_plot = df.set_index(['x', 'y'])

        try:
            vx_2D = df.pivot(index='x', columns='y', values='vx').values
            vy_2D = df.pivot(index='x', columns='y', values='vy').values
            density_2D = df.pivot(index='x', columns='y',
                                  values='density').values
            pressure_2D = df.pivot(
                index='x', columns='y', values='pressure').values
        except:
            # Если названия колонок не заданы, используем индексы
            vx_2D = df.pivot(
                index=df.columns[0], columns=df.columns[1], values=df.columns[2]).values
            vy_2D = df.pivot(
                index=df.columns[0], columns=df.columns[1], values=df.columns[3]).values
            density_2D = df.pivot(
                index=df.columns[0], columns=df.columns[1], values=df.columns[4]).values
            pressure_2D = df.pivot(
                index=df.columns[0], columns=df.columns[1], values=df.columns[5]).values

        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'Time: {t:.2f} s')

        # Тепловая карта скорости по x
        im1 = axs[0, 0].pcolormesh(X, Y, vx_2D, shading='auto', cmap='RdBu_r')
        axs[0, 0].set_title('Velocity X')
        axs[0, 0].set_ylabel('Y Coordinate')
        plt.colorbar(im1, ax=axs[0, 0])

        # Тепловая карта скорости по y
        im2 = axs[0, 1].pcolormesh(X, Y, vy_2D, shading='auto', cmap='RdBu_r')
        axs[0, 1].set_title('Velocity Y')
        plt.colorbar(im2, ax=axs[0, 1])

        # Тепловая карта плотности
        im3 = axs[1, 0].pcolormesh(
            X, Y, density_2D, shading='auto', cmap='viridis')
        axs[1, 0].set_title('Density')
        axs[1, 0].set_xlabel('X Coordinate')
        axs[1, 0].set_ylabel('Y Coordinate')
        plt.colorbar(im3, ax=axs[1, 0])

        # Тепловая карта давления
        im4 = axs[1, 1].pcolormesh(
            X, Y, pressure_2D, shading='auto', cmap='plasma')
        axs[1, 1].set_title('Pressure')
        axs[1, 1].set_xlabel('X Coordinate')
        plt.colorbar(im4, ax=axs[1, 1])

        plt.tight_layout()

        # Сохраняем кадр
        filename = f'frame_2D_{i:03d}.png'
        plt.savefig(filename, dpi=100, bbox_inches='tight')
        plt.close()
        filenames_write.append(filename)

    # Создаем GIF из кадров
    with imageio.get_writer(output_gif, mode='I', duration=0.5) as writer:
        for filename in filenames_write:
            image = imageio.imread(filename)
            writer.append_data(image)

    # Удаляем временные файлы кадров
    for filename in filenames_write:
        os.remove(filename)
