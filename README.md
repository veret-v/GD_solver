
# [AscllGraphic](https://github.com/veret-v/AscllGraphic.git)
---
## Автор

[veret-v](https://github.com/veret-v)

## Стек технологий

[![c++](https://img.shields.io/badge/C++-v13.1.6-red)]()
[![c++](https://img.shields.io/badge/ncurses-v2.6.32-green)]()

## Общий обзор 
---
AscllGraphic это пакет для визуализации 3d моделей, создания сцен с несколькими обьектами, освещением и анимациями. Взаимодействие с программой происходит через терминал linux и json файлы конфигурации. Изображение генерируется в терминале.

## Установка
---
Клонировать репозиторий и перейти в него в командной строке
```terminal
git clone https://github.com/user/ascll_graphic.git
cd ascll_graphic 
```
Собрать бинарный файл
```terminal
mkdir build
cd build
cmake ..
make
```

## Структура пакета
---
```bash
├── CMakeLists.txt
├── CMakePresets.json
├── CMakeUserPresets.json
├── README.md
├── configs
│   ├── config.json
│   └── sphere.json
├── motion
│   ├── circle.json
│   └── line.json
├── objects
│   ├── 1poligon.json
│   ├── 1poligon_points.csv
│   ├── 1poligon_poligons.csv
│   ├── conus.json
│   ├── conus_points.csv
│   ├── conus_poligons.csv
│   ├── object_generators
│   │   ├── conus_generator.cpp
│   │   └── shere_generator.cpp
│   ├── piramid.json
│   ├── piramid_points.csv
│   ├── piramid_poligons.csv
│   ├── plane.json
│   ├── plane_points.csv
│   ├── plane_poligons.csv
│   ├── sphere.json
│   ├── sphere_points.csv
│   ├── sphere_poligons.csv
│   └── tor.json
├── src
│   ├── linear_field_operations
│   │   ├── CMakeLists.txt
│   │   ├── linear_transform.cpp
│   │   ├── linear_transform.h
│   │   ├── point.cpp
│   │   ├── point.h
│   │   ├── vector.cpp
│   │   └── vector.h
│   ├── render
│   │   ├── CMakeLists.txt
│   │   ├── frame.cpp
│   │   ├── frame.h
│   │   ├── render.cpp
│   │   ├── render.h
│   │   ├── vector_2d.cpp
│   │   └── vector_2d.h
│   ├── run_config.cpp
│   └── scene
│       ├── CMakeLists.txt
│       ├── light_source.cpp
│       ├── light_source.h
│       ├── object.cpp
│       ├── object.h
│       ├── point_of_view.cpp
│       ├── point_of_view.h
│       ├── poligon.cpp
│       └── poligon.h
├── vcpkg-configuration.json
└── vcpkg.json

```
* ```code``` содержит исходники проекта
* ```configs``` конфигурационные файлы сцен
* ```motion``` траектории движения камеры и объектов
* ```objects``` объекты которые можно расположить на сцене

## Интерфейс взаимодействия
---
Создание сцен и анимаций происходит через  заполнение конфигурационного файла *config.json*. После чего вводится команда:
```
.\run_scene args
``` 
Ключи ```args``` реализуют дополнительное взаимодействие с выводимым контентом.
Возможные ключи:

* ```config_name``` в случае, если *config.json* не один необходимо ввести этот ключ и название нужного файла, иначе будет выбран первый *.json* файл в папке */config* 
### Структура *config.json*
Содержание поля ```"pov"```  генерируемого контента определяется полем ```"type"``` и может содержать:
* POV
* Движение POV по заданной траектории
```json
{
    "objects": [
        {
            "pos": {
                "x": "double value",
                "y": "double value",
                "z": "double value"
            },
            "scale": "double value",
            "template": "../objects/object.json"
        },
        {
            "pos": {
                "x": "double value",
                "y": "double value",
                "z": "double value"
            },
            "scale": "double value",
            "template": "../objects/object1.json",
            "material": {
                "spec_reflect": 0,
                "diffuse_reflect": 0,
                "specular_exponent": 0,
                "specularity": 0.5,
                "attenuation": [1.0, 0.1, 0.0]
            }
        },
        ...
    ],
    "light_sources": [
        {
            "pos": {
                "x": "double value",
                "y": "double value",
                "z": "double value"
            },
            "intensity": "intensity value int"
        },
        {
            "pos": {
                "x": "double value",
                "y": "double value",
                "z": "double value"
            },
            "intensity": "intensity value int"
        }
        ...
    ],
    "pov": {
        "position": {
                "x": -12.0,
                "y": 0,
                "z": 0
        },
        "normal": {
                "x": 1,
                "y": 0,
                "z": 0
        },
        "tang1": {
                "x": 0,
                "y": 1,
                "z": 0
        },
        "tang2": {
                "x": 0,
                "y": 0,
                "z": 1
        }
    }
}
``` 

### Задание собственных 3D объектов
Создание пользовательских объктов также реализуется путем задания значений в *.json* файле.
Каждый объект состоит из набора полигонов, кодируемых тремя точками в системе отсчета связанной с нулевой системой отсчета .
Примеры лежат в /objects.