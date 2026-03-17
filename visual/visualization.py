import argparse
import glob
import os
import re
import shutil
from pathlib import Path

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METHOD_NAMES = {
    0: "Godunov",
    1: "Kolgan",
    2: "Rodionov",
    3: "HLL",
    4: "HLLC",
    5: "Rusanov",
    6: "Osher",
    7: "Roe",
    8: "FLIC",
    9: "Mader2DE",
}


def load_profile_config(input_ini_path: Path):
    axis = "x"
    profile_index = -1
    gamma = 1.4
    equation_type = -1

    if input_ini_path.exists():
        in_system = False
        with input_ini_path.open("r", encoding="utf-8", errors="ignore") as f:
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
                key = key.lower()
                if key == "analytic_axis":
                    axis = value.lower()
                elif key == "analytic_profile_index":
                    try:
                        profile_index = int(value)
                    except ValueError:
                        pass
                elif key == "g":
                    try:
                        gamma = float(value)
                    except ValueError:
                        pass
                elif key == "equation_type":
                    try:
                        equation_type = int(value)
                    except ValueError:
                        pass

    if axis not in ("x", "y"):
        axis = "x"

    method_tag = f"eq={equation_type} ({METHOD_NAMES.get(equation_type, 'Unknown')})"
    return axis, profile_index, gamma, method_tag


def clamp(value, lo, hi):
    return max(lo, min(value, hi))


def extract_time(path: Path):
    m = re.search(r"time_([\d.]+)\.csv", path.name)
    return float(m.group(1)) if m else None


def select_1d_slice(df, axis, profile_index):
    if "is_solid" in df.columns:
        df = df[df["is_solid"] < 0.5].copy()

    if axis == "x":
        y_vals = np.sort(df["y"].unique())
        idx = len(y_vals) // 4 if profile_index < 0 else clamp(profile_index, 0, len(y_vals) - 1)
        fixed_value = y_vals[idx]
        slice_df = df[np.isclose(df["y"], fixed_value)].copy()
        coord = "x"
        fixed_name = "y"
    else:
        x_vals = np.sort(df["x"].unique())
        idx = len(x_vals) // 4 if profile_index < 0 else clamp(profile_index, 0, len(x_vals) - 1)
        fixed_value = x_vals[idx]
        slice_df = df[np.isclose(df["x"], fixed_value)].copy()
        coord = "y"
        fixed_name = "x"

    # On duplicate rows, average numeric columns by profile coordinate.
    numeric_cols = [c for c in slice_df.columns if pd.api.types.is_numeric_dtype(slice_df[c])]
    slice_df = slice_df[numeric_cols].groupby(coord, as_index=False).mean().sort_values(coord)
    return slice_df, coord, fixed_name, fixed_value


def create_animation_1d(step_files, axis, profile_index, gamma, method_tag, output_gif, frame_dir, duration, use_markers):
    output_gif.parent.mkdir(parents=True, exist_ok=True)
    frame_dir.mkdir(parents=True, exist_ok=True)

    times = []
    for p in step_files:
        t = extract_time(p)
        if t is not None:
            times.append((t, p))
    times.sort(key=lambda x: x[0])

    frame_files = []
    style_num = "o-" if use_markers else "-"

    for idx, (t, path) in enumerate(times):
        df = pd.read_csv(path)
        required = ["x", "y", "u", "rho", "p"]
        if not all(col in df.columns for col in required):
            print(f"Пропуск t={t}: не хватает колонок для 1D профиля")
            continue
        has_exact = all(col in df.columns for col in ["u_exact", "rho_exact", "p_exact"])
        if has_exact:
            has_exact = np.isfinite(df[["u_exact", "rho_exact", "p_exact"]].to_numpy(dtype=float)).any()

        sl, coord, fixed_name, fixed_value = select_1d_slice(df, axis, profile_index)
        if sl.empty:
            print(f"Пропуск t={t}: срез пустой")
            continue

        coord_vals = sl[coord].to_numpy()
        u_num = sl["u"].to_numpy()
        rho_num = sl["rho"].to_numpy()
        p_num = sl["p"].to_numpy()
        e_num = p_num / ((gamma - 1.0) * rho_num)
        if has_exact:
            u_ex = sl["u_exact"].to_numpy()
            rho_ex = sl["rho_exact"].to_numpy()
            p_ex = sl["p_exact"].to_numpy()
            e_ex = p_ex / ((gamma - 1.0) * rho_ex)

        fig, axs = plt.subplots(4, 1, figsize=(12, 12))
        fig.suptitle(f"Time: {t:.6f} s | {method_tag} | slice {coord} at {fixed_name}={fixed_value:.6f}")

        axs[0].plot(coord_vals, u_num, style_num, label="u numerical", linewidth=1.5)
        if has_exact:
            axs[0].plot(coord_vals, u_ex, "-", label="u exact", linewidth=2)
        axs[0].set_ylabel("Velocity")
        axs[0].grid(True, alpha=0.3)
        axs[0].legend()

        axs[1].plot(coord_vals, rho_num, style_num, label="rho numerical", linewidth=1.5)
        if has_exact:
            axs[1].plot(coord_vals, rho_ex, "-", label="rho exact", linewidth=2)
        axs[1].set_ylabel("Density")
        axs[1].grid(True, alpha=0.3)
        axs[1].legend()

        axs[2].plot(coord_vals, p_num, style_num, label="p numerical", linewidth=1.5)
        if has_exact:
            axs[2].plot(coord_vals, p_ex, "-", label="p exact", linewidth=2)
        axs[2].set_ylabel("Pressure")
        axs[2].grid(True, alpha=0.3)
        axs[2].legend()

        axs[3].plot(coord_vals, e_num, style_num, label="e numerical", linewidth=1.5)
        if has_exact:
            axs[3].plot(coord_vals, e_ex, "-", label="e exact", linewidth=2)
        axs[3].set_ylabel("Internal energy")
        axs[3].set_xlabel(coord)
        axs[3].grid(True, alpha=0.3)
        axs[3].legend()

        plt.tight_layout()
        frame = frame_dir / f"frame_1d_{idx:04d}.png"
        plt.savefig(frame, dpi=100, bbox_inches="tight")
        plt.close()
        frame_files.append(frame)
        print(f"Создан кадр {idx + 1}/{len(times)} для времени {t:.6f}")

    if not frame_files:
        print("Не создано ни одного кадра!")
        return 1

    print("Создание GIF...")
    with imageio.get_writer(output_gif, mode="I", duration=duration) as writer:
        for frame in frame_files:
            writer.append_data(imageio.imread(frame))

    shutil.rmtree(frame_dir, ignore_errors=True)
    print(f"GIF успешно создан: {output_gif}")
    return 0


def main():
    project_root = Path(__file__).resolve().parents[1]

    ap = argparse.ArgumentParser(description="Build 1D animation.gif from step_*_time_*.csv")
    ap.add_argument("--steps-dir", default=str(project_root / "output" / "steps"))
    ap.add_argument("--output-gif", default=str(project_root / "animation.gif"))
    ap.add_argument("--input-ini", default=str(project_root / "configs" / "input.ini"))
    ap.add_argument("--frame-dir", default=None)
    ap.add_argument("--duration", type=float, default=0.5)
    ap.add_argument("--no-markers", action="store_true", help="Use plain lines for faster plotting")
    args = ap.parse_args()

    steps_dir = Path(args.steps_dir).resolve()
    output_gif = Path(args.output_gif).resolve()
    input_ini = Path(args.input_ini).resolve()
    frame_dir = Path(args.frame_dir).resolve() if args.frame_dir else (output_gif.parent / ".frames_tmp")

    step_files = sorted(Path(p) for p in glob.glob(str(steps_dir / "step_*_time_*.csv")))
    if not step_files:
        print(f"step_*_time_*.csv не найдены в {steps_dir}")
        return 1

    axis, profile_index, gamma, method_tag = load_profile_config(input_ini)
    print(f"Найдено {len(step_files)} файлов")
    print(f"Профиль для animation.gif: axis={axis}, index={profile_index}, {method_tag}")

    return create_animation_1d(
        step_files,
        axis,
        profile_index,
        gamma,
        method_tag,
        output_gif,
        frame_dir,
        args.duration,
        use_markers=(not args.no_markers),
    )


if __name__ == "__main__":
    raise SystemExit(main())
