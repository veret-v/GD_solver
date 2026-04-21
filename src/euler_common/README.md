## Euler Common

Общий код для эйлеровых задач группы.

В этой папке лежат ключевые части для заданий 2 и 3:
- `solver.cpp` — HLLC, MUSCL/minmod-реконструкция, 2D CFL, 2D state/flux
- `solver.h` — объявления функций
- `grid.h` — индексы переменных и общие константы

Маркеры для проверки:
- `CHECK: HLLC_SOLVER`
- `CHECK: RECONSTRUCTION`
- `CHECK: STATE_VECTOR_2D`
- `CHECK: CFL_2D`
