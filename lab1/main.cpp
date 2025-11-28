#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <random>

using namespace std;

// Структура для хранения данных
struct DataPoint {
    double t;  // время
    double v;  // напряжение
    double c;  // ток
};

// Класс для линейной интерполяции
class LinearInterpolator {
private:
    vector<double> x_data;
    vector<double> y_data;
    
public:
    LinearInterpolator(const vector<double>& x, const vector<double>& y) 
        : x_data(x), y_data(y) {}
    
    double operator()(double x) const {
        // Если x за пределами данных, возвращаем граничные значения
        if (x <= x_data.front()) return y_data.front();
        if (x >= x_data.back()) return y_data.back();
        
        // Находим индекс для интерполяции
        auto it = lower_bound(x_data.begin(), x_data.end(), x);
        size_t idx = distance(x_data.begin(), it);
        
        if (idx == 0) return y_data[0];
        
        double x0 = x_data[idx - 1];
        double x1 = x_data[idx];
        double y0 = y_data[idx - 1];
        double y1 = y_data[idx];
        
        // Линейная интерполяция
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }
};

struct Model {
    vector<double> params;
    
    int steps;
    double dt;
    LinearInterpolator V_func;

    vector<double> t_grid;
    vector<double> V_grid;
    vector<double> I_grid;

    Model(vector<double>& t_exp, vector<double>& V_exp, vector<double>& I_exp, int steps) 
                : steps(steps), V_func(t_exp, V_exp), params(7) {
        LinearInterpolator I_func(t_exp, I_exp);
        double t0 = t_exp.front();
        double t_end = t_exp.back();
        dt = (t_end - t0) / steps;

        t_grid = get_grid_T(t0, t_end, steps);
        V_grid = get_grid(t_grid, V_func);
        I_grid = get_grid(t_grid, I_func);

        // Параметры новой модели
        double Ron = 49.6;      // Ом
        double Roff = 2170.81;    // Ом
        double mu_v = 1853.54;     // коэффициент легирующей подвижности
        double Vp = 0.1143;        // В
        double Vn = -2.55;       // В
        double D = 290.384;       // нм
        double x0 = 0.78;         // начальное состояние

        params[0] = Ron;
        params[1] = Roff;
        params[2] = mu_v;
        params[3] = Vp;
        params[4] = Vn;
        params[5] = D;
        params[6] = x0;
    }

    // Построение сетки времени
    vector<double> get_grid_T(double t0, double t_end, double steps) {
        vector<double> t_grid;

        double dt = (t_end - t0) / steps;
        for (double t = t0; t <= t_end; t += dt) {
            t_grid.push_back(t);
        }

        return t_grid;
    }

    // Построение сетки по интерполяционной функции
    vector<double> get_grid(const vector<double>& t_grid, const LinearInterpolator& func) {
        vector<double> val_grid;

        for (int i = 0; i < t_grid.size(); i++) {
            double t = t_grid[i];
            double val = func(t);
            val_grid.push_back(val);
        }
        
        return val_grid;
    }

    // Ток через мемристор
    double I_model(double x, double V, const vector<double> &params) {
        double Ron = params[0], Roff = params[1];
        // Сопротивление мемристора
        double R_mem = Ron * x + Roff * (1 - x);
        return V / R_mem;
    }

    // Функция для вычисления производной dx/dt
    double dxdt(double x, double V, const vector<double> &params) {
        double I_mem = I_model(x, V, params);
        double Ron = params[0], mu_v = params[2], Vp = params[3], Vn = params[4], D = params[5];
        
        // Вычисляем производную в зависимости от напряжения
        if (V >= Vp) {
            return (mu_v * Vp / (D * D)) * exp(Ron * I_mem / Vp);
        } else if (V <= Vn) {
            return (mu_v * Vn / (D * D)) * exp(Ron * I_mem / Vn);
        } else {
            return mu_v * Ron * I_mem / (D * D);
        }
    }

    double ode(double t, double x, const vector<double> &params) {
        double V = V_func(t);
        return dxdt(x, V, params);
    }

    // Метод Рунге-Кутта 4 порядка для решения ОДУ
    vector<double> solve_x(const vector<double> &params) {
        vector<double> x_sol(steps + 1);
        x_sol[0] = params[6];
        
        double t = t_grid[0];
        for (int i = 1; i <= steps; i++) {
            double k1 = dt * ode(t, x_sol[i-1], params);
            double k2 = dt * ode(t + dt/2, x_sol[i-1] + k1/2, params);
            double k3 = dt * ode(t + dt/2, x_sol[i-1] + k2/2, params);
            double k4 = dt * ode(t + dt, x_sol[i-1] + k3, params);
            
            x_sol[i] = x_sol[i-1] + (k1 + 2*k2 + 2*k3 + k4) / 6;
            
            // Ограничиваем x в диапазоне [0, 1]
            x_sol[i] = max(0.0, min(1.0, x_sol[i]));
            
            t += dt;
        }

        return x_sol;
    }

    vector<double> get_I_solution(const vector<double> &params) {
        vector<double> x_sol = solve_x(params);
        vector<double> I_sol(steps + 1);

        for (int i = 0; i <= steps; i++) {
            double x = x_sol[i];
            double V = V_grid[i];
            I_sol[i] = I_model(x, V, params);
        }
        return I_sol;
    }

    // Целевая функция
    double objective_function(const vector<double> &params) {
        vector<double> I_sol = get_I_solution(params);
        double sum = 0.0;
        for (int i = 0; i < I_sol.size(); i++) {
            double r = I_sol[i] - I_grid[i];
            sum += r*r;
        }

        return sum;
    }

    // // Логорифмическая
    // double objective_function(const vector<double> &params) {
    //     vector<double> I_sol = get_I_solution(params);
    //     double sum = 0.0;
    //     int count = 0;
        
    //     for (int i = 0; i < I_sol.size(); i++) {
    //         if (I_sol[i] > 1e-12 && I_grid[i] > 1e-12) {
    //             double log_error = log(fabs(I_sol[i])) - log(fabs(I_grid[i]));
    //             sum += log_error * log_error;
    //             count++;
    //         }
    //     }
    //     return count > 0 ? sum / count : 1e12;
    // }
};

// Функция для чтения CSV файла
vector<DataPoint> readCSV(const string& filename) {
    vector<DataPoint> data;
    ifstream file(filename);
    string line;
    
    if (!file.is_open()) {
        cerr << "Ошибка открытия файла: " << filename << endl;
        return data;
    }
    
    // Пропускаем заголовок 
    getline(file, line);
    
    while (getline(file, line)) {
        stringstream ss(line);
        DataPoint point;
        string token;
        
        // Читаем время
        if (!getline(ss, token, '\t')) continue;
        point.t = stod(token);
        
        // Читаем напряжение
        if (!getline(ss, token, '\t')) continue;
        point.v = stod(token);
        
        // Читаем ток
        if (!getline(ss, token, '\t')) continue;
        point.c = stod(token);

        
        data.push_back(point);
    }
    
    return data;
}

// Функция для систематической генерации начальных параметров
vector<vector<double>> generate_parameter_samples(const vector<double>& base_params) {
    vector<vector<double>> samples;
    int N = base_params.size();
    
    // Создаем сетку значений для каждого параметра
    vector<vector<double>> param_values(N);
    
    for (int i = 0; i < N; i++) {
        double base_value = base_params[i];
        vector<double> values;
        
        if (i == 6) { // x0 - особый случай [0, 1]
            values = {0.1, 0.2, 0.5, 0.7, 0.9}; // фиксированные значения для x0
        } else {
            // Систематические значения
            values = {base_value - base_value * 0.1, base_value - base_value * 0.01, base_value, base_value + base_value * 0.01, base_value + base_value * 0.1};
        }
        param_values[i] = values;
    }
    
    // Генерируем все возможные комбинации
    vector<int> indices(N, 0);
    
    while (true) {
        vector<double> sample(N);
        for (int i = 0; i < N; i++) {
            sample[i] = param_values[i][indices[i]];
        }
        samples.push_back(sample);
        
        // Увеличиваем индексы
        int j = N - 1;
        while (j >= 0) {
            indices[j]++;
            if (indices[j] < param_values[j].size()) {
                break;
            }
            indices[j] = 0;
            j--;
        }
        if (j < 0) break;
    }
    
    return samples;
}

// Функция для создания выборки начальных параметров и выбора лучшей
vector<double> find_best_initial_parameters(Model& model) {
    vector<double> best_params = model.params;
    double best_value = model.objective_function(best_params);
    
    cout << "Систематический поиск лучших начальных параметров..." << endl;
    cout << "Исходные параметры: J = " << best_value << endl;
    
    // Генерируем систематическую выборку
    auto samples = generate_parameter_samples(model.params);
    
    cout << "Всего вариантов для проверки: " << samples.size() << endl;
    
    for (size_t i = 0; i < samples.size(); i++) {
        vector<double> candidate_params = samples[i];
        double candidate_value = model.objective_function(candidate_params);
        
        cout << "Вариант " << i + 1 << "/" << samples.size() << ": J = " << candidate_value << endl;
        
        if (candidate_value < best_value) {
            best_value = candidate_value;
            best_params = candidate_params;
            cout << ">>> НАЙДЕН УЛУЧШЕННЫЙ ВАРИАНТ! J = " << best_value << endl;
        }
    }
    
    cout << "Лучшие начальные параметры найдены с J = " << best_value << endl;
   cout << "Ron=" << best_params[0] << ", Roff=" << best_params[1] 
         << ", mu_v=" << best_params[2] << ", Vp=" << best_params[3] 
         << ", Vn=" << best_params[4] << ", D=" << best_params[5] 
         << ", x0=" << best_params[6] << endl << endl;
    
    return best_params;
}

// Исследующий поиск по координатным направлениям (правильный)
void exploratory_search(Model& model, const vector<double>& step_sizes, 
                       vector<double>& current_params, double& current_value) {
    int size_p = model.params.size();
    
    for (int i = 0; i < size_p; i++) {
        // Пробуем шаг в положительном направлении
        vector<double> params_positive = current_params;
        params_positive[i] += step_sizes[i];
        double f_positive = model.objective_function(params_positive);
        
        if (f_positive < current_value) {
            current_params = params_positive;
            current_value = f_positive;
            continue;  // Успех - переходим к следующей координате
        }
        
        // Пробуем шаг в отрицательном направлении
        vector<double> params_negative = current_params;
        params_negative[i] -= step_sizes[i];
        double f_negative = model.objective_function(params_negative);
        
        if (f_negative < current_value) {
            current_params = params_negative;
            current_value = f_negative;
        }
        // Если оба шага неудачны - оставляем параметр без изменений
    }
}

// Основной метод оптимизации Хука-Дживса
void optimize(Model& model, double eps = 1e-6, double init_step = 0.1,
              double accel = 2.0, double reduction = 0.5, int max_iter = 1000) {
    
    int size_p = model.params.size();
    
    // Начальная точка
    vector<double> x_current = model.params;
    //x_current = find_best_initial_parameters(model);
    double f_current = model.objective_function(x_current);
    
    // Базовая точка
    vector<double> x_base = x_current;
    double f_base = f_current;
    
    // Вектор шагов (процент от текущих значений параметров)
    vector<double> step_sizes(size_p);
    for(int i = 0; i < size_p; i++) {
        step_sizes[i] = x_current[i] * 0.1;
    }

    
    int iteration = 0;
    
    while (iteration < max_iter) {
        // Сохраняем точку ДО исследующего поиска
        vector<double> x_before_search = x_current;
        double f_before_search = f_current;
        
        // Исследующий поиск
        exploratory_search(model, step_sizes, x_current, f_current);
        
        // Проверяем, был ли исследующий поиск успешным
        bool search_successful = (f_current < f_before_search - eps);
        
        if (search_successful) {
            // Успешный исследующий поиск - поиск по образцу
            vector<double> x_pattern(size_p);
            for (int i = 0; i < size_p; i++) {
                x_pattern[i] = x_current[i] + accel * (x_current[i] - x_base[i]);
            }
            
            double f_pattern = model.objective_function(x_pattern);
            
            if (f_pattern < f_current) {
                // Поиск по образцу успешен
                x_base = x_current;  // Новая базовая точка
                f_base = f_current;
                x_current = x_pattern;
                f_current = f_pattern;
            } else {
                // Поиск по образцу неудачен
                x_base = x_current;  // Новая базовая точка
                f_base = f_current;
            }
        } else {
            // Исследующий поиск неудачен - уменьшаем шаги
            bool all_steps_small = true;
            for (int i = 0; i < size_p; i++) {
                step_sizes[i] *= reduction;
                if (step_sizes[i] > eps) {
                    all_steps_small = false;
                }
            }
            
            // Если все шаги стали очень маленькими - завершаем
            if (all_steps_small) {
                cout << "Оптимизация завершена: шаги стали слишком малыми" << endl;
                break;
            }
            
            // Возвращаемся к базовой точке для следующего поиска
            x_current = x_base;
            f_current = f_base;
        }
        
        // Обновляем параметры модели
        model.params = x_current;
        
        // Вывод информации
        if (iteration % 10 == 0) {
            cout << "Итерация " << iteration << ": f = " << f_current << endl;
            cout << "Шаги: ";
            for (double step : step_sizes) cout << step << " ";
            cout << endl;
        }
        
        iteration++;
    }
    
    cout << "Финальный результат:" << endl;
    cout << "f = " << f_current << " после " << iteration << " итераций" << endl;
   cout << "Ron=" << model.params[0] << ", Roff=" << model.params[1] 
         << ", mu_v=" << model.params[2] << ", Vp=" << model.params[3] 
         << ", Vn=" << model.params[4] << ", D=" << model.params[5] 
         << ", x0=" << model.params[6] << endl << endl;
}

int main() {
    ///////////////////////////////////////////////////////////////// Чтение данных ////////////////////////////////////////////////////////////////////////////

    auto data = readCSV("data/LiNbO3.csv");
    
    if (data.empty()) {
        cerr << "Нет данных для обработки" << endl;
        return 1;
    }
    
    cout << "Прочитано " << data.size() << " точек данных" << endl;
    cout << "Первые 5 строк:" << endl;
    for (int i = 0; i < min(5, (int)data.size()); i++) {
        cout << "t=" << data[i].t << ", v=" << data[i].v << ", c=" << data[i].c << endl;
    }
    cout << endl;

    // Извлекаем времена и напряжения
    vector<double> t_exp, V_exp, I_exp;
    for (const auto& point : data) {
        t_exp.push_back(point.t);
        V_exp.push_back(point.v);
        I_exp.push_back(point.c);
    }

    ////////////////////////////////////////////////////////////////////// Расчет ////////////////////////////////////////////////////////////////////////////
    int steps = 5000;
    Model model(t_exp, V_exp, I_exp, steps);

    cout << "Начинаем поиск параметров модели..." << endl;
    cout << "Ошибка: " << model.objective_function(model.params) << endl;

    optimize(model);

    // Сохраняем результаты
    ofstream out_file("memristor_results.csv");

    out_file << "t,V,I\n";
    cout << "Ron=" << model.params[0] << ", Roff=" << model.params[1] 
         << ", mu_v=" << model.params[2] << ", Vp=" << model.params[3] 
         << ", Vn=" << model.params[4] << ", D=" << model.params[5] 
         << ", x0=" << model.params[6] << endl << endl;
         
    vector<double> I_sol = model.get_I_solution(model.params);
    cout << "Ошибка: " << model.objective_function(model.params) << endl;
    
    for (size_t i = 0; i <= steps; i++) {
        out_file << model.t_grid[i] << "," 
            << model.V_grid[i] << ","
            << I_sol[i] << "\n";
    }
    
    out_file.close();
    
    cout << "Результаты сохранены в memristor_results.csv" << endl;
    cout << "Количество точек расчета: " << model.t_grid.size() << endl;

    return 0;
}