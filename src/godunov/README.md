## Godunov

Конфиги и привязка к результатам для заданий Godunov 1D/2D.

Файлы:
- `sod_hllc_1d.ini` — запуск 1D Sod для `results/02_sod_godunov/`
- `sod_hllc_2d.ini` — 2D экструзия Sod для `results/03_godunov_2d/extruded_1d.png`
- `conservation_2d.ini` — отдельный консервативный 2D тест для `conservation.txt`
- `vortex_N40.ini`, `vortex_N80.ini`, `vortex_N160.ini` — сетки для `convergence.png`

Численный общий код для этих конфигов лежит в `../euler_common/`.
