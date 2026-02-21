import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import glob
import imageio
import os
import re

def create_animation_1D(filenames, output_gif='animation.gif'):
    # Извлекаем время из названий файлов и создаем словарь данных
    data = {}
    for filename in filenames:
        # Используем регулярное выражение для извлечения времени
        match = re.search(r'time_([\d.]+)\.csv', filename)
        if match:
            time = float(match.group(1))
            df = pd.read_csv(filename)
            data[time] = df
    
    # Сортируем временные точки
    times = sorted(data.keys())
    
    filenames_write = []
    
    # Создаем кадры для каждого момента времени
    for i, t in enumerate(times):
        df = data[t]
        
        # Проверяем наличие необходимых колонок
        if not all(col in df.columns for col in ['x', 'u', 'rho', 'p']):
            print(f"Предупреждение: в файле для времени {t} отсутствуют некоторые колонки")
            print(f"Доступные колонки: {df.columns.tolist()}")
            continue
        
        # Извлекаем данные
        x_coord = df['x']
        velocity = df['u']
        velocity_exact = df['u_exact']
        density = df['rho']
        density_exact = df['rho_exact']
        pressure = df['p']
        pressure_exact = df['p_exact']
        
        fig, axs = plt.subplots(4, 1, figsize=(12, 12))
        fig.suptitle(f'Time: {t:.6f} s')
        
        # График скорости от координаты
        axs[0].plot(x_coord, velocity, 'o-', label='analitycal_velocity', linewidth=2)
        axs[0].plot(x_coord, velocity_exact, 'b-', label='numerical_velocity', linewidth=2)
        axs[0].set_ylabel('Velocity')
        axs[0].grid(True, alpha=0.3)
        axs[0].legend()
        
        # График плотности от координаты
        axs[1].plot(x_coord, density, 'o-', label='Density_numerical', linewidth=2)
        axs[1].plot(x_coord, density_exact, 'b-', label='Density_analitycal', linewidth=2)
        axs[1].set_ylabel('Density')
        axs[1].grid(True, alpha=0.3)
        axs[1].legend()
        
        # График давления от координаты
        axs[2].plot(x_coord, pressure, 'o-', label='Pressure_numerical', linewidth=2)
        axs[2].plot(x_coord, pressure_exact, 'b-', label='Pressure_analitycal', linewidth=2)
        axs[2].set_ylabel('Pressure')
        axs[2].set_xlabel('Coordinate')
        axs[2].grid(True, alpha=0.3)
        axs[2].legend()
        
        internal_energy = (1 / 0.4) * pressure / density
        internal_energy_exact = (1 / 0.4) * pressure_exact / density_exact
        axs[3].plot(x_coord, internal_energy, 'o-', label='internal_energy', linewidth=2)
        axs[3].plot(x_coord, internal_energy_exact, 'b-', label='internal_energy_exact', linewidth=2)
        axs[3].set_ylabel('Internal energy')
        axs[3].set_xlabel('Coordinate')
        axs[3].grid(True, alpha=0.3)
        axs[3].legend()
        
        plt.tight_layout()
        
        # Сохраняем кадр
        filename = f'frame_{i:03d}.png'
        plt.savefig(filename, dpi=100, bbox_inches='tight')
        plt.close()
        filenames_write.append(filename)
        
        print(f"Создан кадр {i+1}/{len(times)} для времени {t}")
    
    # Создаем GIF из кадров
    if filenames_write:
        print("Создание GIF...")
        with imageio.get_writer(output_gif, mode='I', duration=0.5) as writer:
            for filename in filenames_write:
                image = imageio.imread(filename)
                writer.append_data(image)
        
        # Удаляем временные файлы кадров
        for filename in filenames_write:
            os.remove(filename)
        
        print(f"GIF успешно создан: {output_gif}")
    else:
        print("Не создано ни одного кадра!")

# Основной код
if __name__ == "__main__":
    # Исправленный путь - используем сырую строку или двойные обратные слеши
    path = r'C:\solvver\GD_solver\output\steps'  # Сырая строка
    # Или: path = 'C:\\solvver\\GD_solver\\output\\steps'  # Двойные обратные слеши
    
    # Проверяем существование пути
    if not os.path.exists(path):
        print(f"Ошибка: путь {path} не существует!")
        exit()
    
    # Ищем файлы
    filenames = glob.glob(os.path.join(path, "*.csv"))
    
    if not filenames:
        print(f"CSV файлы не найдены в папке: {path}")
        # Попробуем найти файлы в текущей директории
        filenames = glob.glob("step_*_time_*.csv")
        if filenames:
            print("Найдены файлы в текущей директории")
    
    if not filenames:
        print("Файлы не найдены!")
        exit()
    
    print(f"Найдено {len(filenames)} файлов")
    
    # Показываем структуру первого файла для отладки
    if filenames:
        sample_df = pd.read_csv(filenames[0])
        print(f"Структура первого файла ({filenames[0]}):")
        print(f"Колонки: {sample_df.columns.tolist()}")
        print(f"Размер: {sample_df.shape}")
    
    create_animation_1D(filenames)