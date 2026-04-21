# Makefile для GodunovSolver (C++17, оптимизации)
# Компилятор: g++
CXX = g++

# Стандарт C++ и оптимизации
CXXFLAGS = -std=c++17 -O3 

# Предупреждения (полезно при разработке, можно убрать)
CXXFLAGS += -Wall -Wextra

# Директории для поиска заголовочных файлов
INCLUDES = -I. -I./src

# Исходные файлы (все .cpp из src)
SRCS = src/main.cpp \
       src/parser.cpp \
       src/point.cpp \
       src/grid.cpp \
       src/solver.cpp \
       src/script.cpp

# Объектные файлы (заменяем .cpp на .o, кладём в текущую папку)
OBJS = $(SRCS:src/%.cpp=%.o)

# Имя исполняемого файла
TARGET = solver

# Цель по умолчанию: сборка программы
all: $(TARGET)

# Линковка
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Компиляция .cpp в .o с автоматической генерацией зависимостей
%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

# Подключение сгенерированных зависимостей (файлов .d)
-include $(OBJS:.o=.d)

# Очистка
clean:
	rm -f $(TARGET) $(OBJS) $(OBJS:.o=.d)

# Пересборка с нуля
distclean: clean
	rm -f *~

# Защита от случайного удаления
.PHONY: all clean distclean