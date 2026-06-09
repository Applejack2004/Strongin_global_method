import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import interp1d

# Пути к файлам
data_file_path = "graphics.txt"
stats_file_path = "stats.txt"

# Чтение данных
try:
    data = np.loadtxt(data_file_path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
except Exception as e:
    print(f"Error reading graphics.txt: {e}")
    data = np.empty((0, 2))

try:
    with open(stats_file_path, 'r') as f:
        line = f.readline().split()
        if len(line) >= 3:
            epsilon, r, numTests = map(float, line[:3])
        else:
            raise ValueError("Not enough arguments in stats.txt")
except Exception as e:
    print(f"Error reading stats file: {e}")
    # Значения по умолчанию
    epsilon, r, numTests = (0.01, 2.0, len(data) if data.size > 0 else 0)

# Обработка данных
if data.size > 0:
    # data[:, 1] это количество итераций
    iterations = data[:, 1]
    
    # Общее число задач берем из stats.txt
    total_problems = int(numTests) if numTests > 0 else len(iterations)
    
    # Сортируем итерации по возрастанию
    iterations_sorted = np.sort(iterations)
    
    # Формируем ось Y: накопительный процент решенных задач
    # np.arange(1, ...) создает массив [1, 2, 3, ... N]
    # Делим на общее число задач и умножаем на 100
    solved_percentage = (np.arange(1, len(iterations_sorted) + 1) / total_problems) * 100.0
else:
    iterations_sorted = np.array([])
    solved_percentage = np.array([])
    total_problems = 1

# Среднее число итераций
avg_iterations = np.mean(iterations_sorted) if len(iterations_sorted) > 0 else 0

# Функция интерполяции (для красивой линии)
def get_interpolation(x, y):
    if len(x) < 2:
        return x, y
    # Интерполяция для более плавного графика, если точек много
    f = interp1d(x, y, kind='linear', bounds_error=False, fill_value=(y[0], y[-1]))
    x_new = np.linspace(x.min(), x.max(), 500)
    y_new = f(x_new)
    return x_new, y_new

# Построение графика
def plot_graph():
    plt.figure(figsize=(9, 7))
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    if iterations_sorted.size > 0:
        # Рисуем точки и соединяем их линией
        plt.plot(iterations_sorted, solved_percentage, 
                 marker='o', markersize=4, linestyle='-', 
                 label='Experimental Data', color='royalblue')
        
    plt.title(f"Operational Characteristics\n(epsilon={epsilon}, r={r}, numTests={int(numTests)})")
    plt.ylabel("Percentage of Solved Problems (%)")
    plt.xlabel("Number of Iterations")
    
    # !!! ЖЕСТКАЯ ФИКСАЦИЯ ОСИ Y ОТ 0 ДО 105 (чтобы 100 влезало красиво) !!!
    plt.ylim(0, 105)
    
    # Добавляем легенду и статистику
    plt.legend(loc='lower right')
    
    if avg_iterations > 0:
        info_text = f'Avg iterations: {avg_iterations:.2f}\nSolved: {len(iterations_sorted)}/{int(numTests)}'
        plt.text(0.02, 0.95, info_text, transform=plt.gca().transAxes, 
                 fontsize=10, verticalalignment='top', 
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if iterations_sorted.size == 0:
        print("No data available to plot.")
    else:
        plot_graph()