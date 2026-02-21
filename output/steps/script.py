import os
import re
import glob
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import imageio
from matplotlib.colors import Normalize

def extract_time(filename):
    """Извлекает значение времени из имени файла."""
    match = re.search(r'time_([0-9]+\.?[0-9]*)', filename)
    if match:
        return float(match.group(1))
    return 0.0

def prepare_grid(df, var_name):
    """Преобразует таблицу с координатами в 2D матрицу."""
    x_vals = np.sort(df['x'].unique())
    y_vals = np.sort(df['y'].unique())
    nx, ny = len(x_vals), len(y_vals)
    mat = np.full((nx, ny), np.nan)
    for i, x in enumerate(x_vals):
        for j, y in enumerate(y_vals):
            val = df[(df['x'] == x) & (df['y'] == y)][var_name]
            if not val.empty:
                mat[i, j] = val.iloc[0]
    return x_vals, y_vals, mat

def create_frame_2d(df, time, var_name, vmin=None, vmax=None, output_path=None):
    """Создаёт один кадр (2D-поле) и сохраняет его."""
    x, y, Z = prepare_grid(df, var_name)
    
    # Защита от vmin == vmax
    if vmin is not None and vmax is not None and vmin == vmax:
        vmin -= 1e-6
        vmax += 1e-6
    
    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(x, y, Z.T, shading='auto', cmap='viridis',
                         norm=Normalize(vmin=vmin, vmax=vmax))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'{var_name} at time = {time:.4f}')
    fig.colorbar(mesh, ax=ax, label=var_name)
    fig.tight_layout()
    fig.savefig(output_path, dpi=100)
    plt.close(fig)

def create_profile_x_frame(df, time, var_name, y_fixed, vmin, vmax, output_path, exact_col=None):
    """
    Создаёт кадр профиля вдоль оси X при фиксированном y = y_fixed.
    Если exact_col задан и существует в df, добавляется кривая точного решения.
    """
    # Находим ближайшее доступное значение y
    unique_y = df['y'].unique()
    nearest_y = unique_y[np.argmin(np.abs(unique_y - y_fixed))]
    subset = df[np.isclose(df['y'], nearest_y)].sort_values('x')
    
    fig, ax = plt.subplots(figsize=(8, 4))
    # Численное решение
    ax.plot(subset['x'], subset[var_name], 'b-', label='numerical')
    # Точное решение, если есть
    if exact_col and exact_col in df.columns:
        ax.plot(subset['x'], subset[exact_col], 'k--', label='exact')
        ax.legend(loc='best')
    
    ax.set_xlabel('x')
    ax.set_ylabel(var_name)
    ax.set_title(f'{var_name} profile at y = {nearest_y:.3f}, time = {time:.4f}')
    ax.set_ylim(vmin, vmax)
    ax.grid(True, linestyle='--', alpha=0.7)
    fig.tight_layout()
    fig.savefig(output_path, dpi=100)
    plt.close(fig)

def create_profile_y_frame(df, time, var_name, x_fixed, vmin, vmax, output_path, exact_col=None):
    """
    Создаёт кадр профиля вдоль оси Y при фиксированном x = x_fixed.
    Если exact_col задан и существует в df, добавляется кривая точного решения.
    """
    unique_x = df['x'].unique()
    nearest_x = unique_x[np.argmin(np.abs(unique_x - x_fixed))]
    subset = df[np.isclose(df['x'], nearest_x)].sort_values('y')
    
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(subset['y'], subset[var_name], 'r-', label='numerical')
    if exact_col and exact_col in df.columns:
        ax.plot(subset['y'], subset[exact_col], 'k--', label='exact')
        ax.legend(loc='best')
    
    ax.set_xlabel('y')
    ax.set_ylabel(var_name)
    ax.set_title(f'{var_name} profile at x = {nearest_x:.3f}, time = {time:.4f}')
    ax.set_ylim(vmin, vmax)
    ax.grid(True, linestyle='--', alpha=0.7)
    fig.tight_layout()
    fig.savefig(output_path, dpi=100)
    plt.close(fig)

def create_gif(file_pattern, var_name, output_gif='animation.gif', fps=5,
               make_2d=True, make_profile_x=False, make_profile_y=False,
               y_fixed=None, x_fixed=None, exact_suffix='_exact'):
    """
    Основная функция: читает файлы, создаёт кадры (2D и/или профили) и собирает GIF.
    
    Параметры:
        file_pattern   : шаблон для поиска CSV-файлов.
        var_name       : имя переменной для визуализации.
        output_gif     : базовое имя для выходного GIF (для профилей будут добавлены суффиксы).
        fps            : кадров в секунду.
        make_2d        : создавать ли анимацию 2D-поля.
        make_profile_x : создавать ли анимацию профиля по X.
        make_profile_y : создавать ли анимацию профиля по Y.
        y_fixed        : фиксированное значение y для профиля X (если None – берётся середина).
        x_fixed        : фиксированное значение x для профиля Y (если None – берётся середина).
        exact_suffix   : суффикс для имени колонки точного решения (например, '_exact').
                         Итоговая колонка будет var_name + exact_suffix.
    """
    files = glob.glob(file_pattern)
    if not files:
        print(f"Файлы по шаблону '{file_pattern}' не найдены.")
        return
    
    files_with_time = [(f, extract_time(f)) for f in files]
    files_with_time.sort(key=lambda x: x[1])
    
    exact_col = var_name + exact_suffix  # имя колонки с точным решением
    
    # Сбор данных для глобальных min/max (включая точные решения, если они есть)
    all_data = []
    valid_files = []   # (filename, time, DataFrame)
    for f, t in files_with_time:
        try:
            df = pd.read_csv(f, float_precision='high')
            if var_name not in df.columns:
                print(f"Предупреждение: в файле {f} нет колонки '{var_name}', пропускаем.")
                continue
            # Добавляем численные значения
            all_data.append(df[var_name].values)
            # Если есть точное решение, добавляем и его
            if exact_col in df.columns:
                all_data.append(df[exact_col].values)
            valid_files.append((f, t, df))
        except Exception as e:
            print(f"Ошибка при чтении файла {f}: {e}")
            continue
    
    if not all_data:
        print("Нет данных для обработки.")
        return
    
    global_min = np.min(np.concatenate(all_data))
    global_max = np.max(np.concatenate(all_data))
    print(f"Глобальный диапазон (с учётом точных решений): [{global_min:.4f}, {global_max:.4f}]")
    
    # Определяем фиксированные координаты для профилей (если не заданы)
    if make_profile_x and y_fixed is None:
        sample_df = valid_files[0][2]
        y_vals = sample_df['y'].unique()
        y_fixed = np.mean([y_vals.min(), y_vals.max()])
        print(f"y_fixed для профиля X не задан, используем {y_fixed:.3f}")
    
    if make_profile_y and x_fixed is None:
        sample_df = valid_files[0][2]
        x_vals = sample_df['x'].unique()
        x_fixed = np.mean([x_vals.min(), x_vals.max()])
        print(f"x_fixed для профиля Y не задан, используем {x_fixed:.3f}")
    
    # Создание временной папки
    temp_dir = 'frames_temp'
    os.makedirs(temp_dir, exist_ok=True)
    
    # Списки для хранения путей к кадрам каждого типа
    frames_2d = []
    frames_prof_x = []
    frames_prof_y = []
    
    for i, (f, time_val, df) in enumerate(valid_files):
        print(f"Обработка {i+1}/{len(valid_files)}: {os.path.basename(f)}")
        
        if make_2d:
            frame_path = os.path.join(temp_dir, f"frame_2d_{i:04d}.png")
            try:
                create_frame_2d(df, time_val, var_name, global_min, global_max, frame_path)
                frames_2d.append(frame_path)
            except Exception as e:
                print(f"Ошибка при создании 2D-кадра для {f}: {e}")
        
        if make_profile_x:
            frame_path = os.path.join(temp_dir, f"frame_prof_x_{i:04d}.png")
            try:
                create_profile_x_frame(df, time_val, var_name, y_fixed,
                                        global_min, global_max, frame_path,
                                        exact_col=exact_col if exact_col in df.columns else None)
                frames_prof_x.append(frame_path)
            except Exception as e:
                print(f"Ошибка при создании профиля X для {f}: {e}")
        
        if make_profile_y:
            frame_path = os.path.join(temp_dir, f"frame_prof_y_{i:04d}.png")
            try:
                create_profile_y_frame(df, time_val, var_name, x_fixed,
                                        global_min, global_max, frame_path,
                                        exact_col=exact_col if exact_col in df.columns else None)
                frames_prof_y.append(frame_path)
            except Exception as e:
                print(f"Ошибка при создании профиля Y для {f}: {e}")
    
    # Функция для сборки GIF из списка кадров
    def assemble_gif(frame_list, gif_name):
        if not frame_list:
            print(f"Нет кадров для {gif_name}, пропускаем.")
            return
        print(f"Создание GIF: {gif_name} (fps={fps})")
        try:
            with imageio.get_writer(gif_name, mode='I', fps=fps) as writer:
                for frame_path in frame_list:
                    if not os.path.exists(frame_path) or os.path.getsize(frame_path) == 0:
                        print(f"Предупреждение: файл {frame_path} пропущен.")
                        continue
                    image = imageio.imread(frame_path)
                    writer.append_data(image)
            print(f"GIF успешно создан: {gif_name}")
        except Exception as e:
            print(f"Ошибка при создании GIF {gif_name}: {e}")
    
    # Собираем GIF-файлы
    base, ext = os.path.splitext(output_gif)
    if make_2d:
        assemble_gif(frames_2d, output_gif)
    if make_profile_x:
        assemble_gif(frames_prof_x, f"{base}_profile_x{ext}")
    if make_profile_y:
        assemble_gif(frames_prof_y, f"{base}_profile_y{ext}")
    
    # Очистка временной папки
    try:
        shutil.rmtree(temp_dir)
        print(f"Временная папка {temp_dir} удалена.")
    except Exception as e:
        print(f"Не удалось удалить временную папку: {e}")

if __name__ == "__main__":
    # Пример: для 'rho' создаём 2D, профили X и Y с точным решением (колонка 'rho_exact')
    create_gif('step_*.csv', var_name='rho',
               output_gif='rho_evolution.gif', fps=5,
               make_2d=True, make_profile_x=True, make_profile_y=True,
               exact_suffix='_exact')
    
    # Для других переменных можно аналогично:
    create_gif('step_*.csv', var_name='u',
                output_gif='u_evolution.gif', fps=5,
                make_2d=True, make_profile_x=True, make_profile_y=True,
                exact_suffix='_exact')
    
    create_gif('step_*.csv', var_name='v',
                output_gif='v_evolution.gif', fps=5,
                make_2d=True, make_profile_x=True, make_profile_y=True,
                exact_suffix='_exact')
    
    create_gif('step_*.csv', var_name='p',
                output_gif='p_evolution.gif', fps=5,
                make_2d=True, make_profile_x=True, make_profile_y=True,
                exact_suffix='_exact')
    