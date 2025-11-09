import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
DEFAULT_CSV_PATH = '/home/dmitrytorov/deniska/GD_solver/output/final_results.csv'


def plot_final_state(csv_path: str = DEFAULT_CSV_PATH, save_path: str | None = None):
    """
    Строит три графика (скорость, плотность, давление) по данным из итогового CSV.

    Parameters
    ----------
    csv_path : str
        Полный путь к файлу результатов.
    save_path : str | None
        Если указан, сохранит рисунок в файл. Иначе покажет окно matplotlib.
    """
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path)
    expected_cols = ['x', 'rho', 'u', 'p']
    if list(df.columns[:4]) != expected_cols:
        raise ValueError(f"Unexpected columns in CSV: {df.columns.tolist()}. "
                         f"Expected at least {expected_cols}.")

    x = df['x']
    rho = df['rho']
    velocity = df['u']
    pressure = df['p']

    fig, axs = plt.subplots(3, 1, figsize=(4, 4), sharex=True)
    fig.suptitle('Final solution profiles')

    axs[0].plot(x, velocity, color='tab:blue', linewidth=2)
    axs[0].set_ylabel('Velocity')
    axs[0].grid(True, alpha=0.3)

    axs[1].plot(x, rho, color='tab:green', linewidth=2)
    axs[1].set_ylabel('Density')
    axs[1].grid(True, alpha=0.3)

    axs[2].plot(x, pressure, color='tab:red', linewidth=2)
    axs[2].set_ylabel('Pressure')
    axs[2].set_xlabel('Coordinate x')
    axs[2].grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if save_path:
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
    else:
        plt.show()


if __name__ == "__main__":
    plot_final_state()