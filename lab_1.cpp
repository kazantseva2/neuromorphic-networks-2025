#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>

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

    vector<double> x_sol;
    vector<double> I_sol;


    Model(vector<double>& t_exp, vector<double>& V_exp, vector<double>& I_exp, int steps) : steps(steps), V_func(t_exp, V_exp), x_sol(steps + 1), I_sol(steps + 1), params(7) {
        LinearInterpolator I_func(t_exp, I_exp);
        double t0 = t_exp.front();
        double t_end = t_exp.back();
        dt = (t_end - t0) / steps;

        t_grid = get_grid_T(t0, t_end, steps);
        V_grid = get_grid(t_grid, V_func);
        I_grid = get_grid(t_grid, I_func);

        // Параметры мемристора
        // double Ron = 205.0;      // Ом
        // double Roff = 2130.0;    // Ом
        // double mu_v = 1000;     // коэффициент легирующей подвижности
        // double Vp = 0.65;        // В
        // double Vn = -0.87;       // В
        // double D = 620;       // нм
        // double x0 = 0.1;         // начальное состояние

        double Ron = 40.0;      // Ом
        double Roff = 2130.0;    // Ом
        double mu_v = 1034;     // коэффициент легирующей подвижности
        double Vp = 0.40;        // В
        double Vn = -1.25;       // В
        double D = 584;       // нм
        double x0 = 1.1;         // начальное состояние

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
    double I(double x, double V) {
        double Ron = params[0], Roff = params[1];
        // Сопротивление мемристора
        double R_mem = Ron * x + Roff * (1 - x);
        return V / R_mem;
    }

    // Функция для вычисления производной dx/dt
    double dxdt(double x, double V) {
        double I_mem = I(x, V);
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

    double ode(double t, double x) {
        double V = V_func(t);
        return dxdt(x, V);
    }
   
    // Метод Рунге-Кутта 4 порядка для решения ОДУ
    void solve_x() {
        x_sol[0] = params[6];
        
        double t = t_grid[0];
        for (int i = 1; i <= steps; i++) {
            double k1 = dt * ode(t, x_sol[i-1]);
            double k2 = dt * ode(t + dt/2, x_sol[i-1] + k1/2);
            double k3 = dt * ode(t + dt/2, x_sol[i-1] + k2/2);
            double k4 = dt * ode(t + dt, x_sol[i-1] + k3);
            
            x_sol[i] = x_sol[i-1] + (k1 + 2*k2 + 2*k3 + k4) / 6;
            
            // Ограничиваем x в диапазоне [0, 1]
            x_sol[i] = max(0.0, min(1.0, x_sol[i]));
            
            t += dt;
        }
    }

    // Функция для обновления векторов решений 
    void update_solution() {
        solve_x();

        for (int i = 0; i <= steps; i++) {
            double x = x_sol[i];
            double V = V_grid[i];
            I_sol[i] = I(x, V);
        }
    }

    // Целевая функция
    double J() {
        update_solution();

        double sum = 0.0;
        for (int i = 0; i < I_sol.size(); i++) {
            double r = I_sol[i] - I_grid[i];
            // r /= I_grid[i];
            sum += r*r;
        }

        return sum;
    }

    double J(vector<double> params_new) {
        vector<double> params_current = params;
        params = params_new;
        update_solution();

        double sum = 0.0;
        for (int i = 0; i < I_sol.size(); i++) {
            double r = I_sol[i] - I_grid[i];
            // r /= I_grid[i];
            sum += r*r;
        }

        params = params_current;

        return sum;
    }

};

// Исследующий поиск по координатным направлениям
void exploratory_search(Model& model, const vector<double>& step_sizes, 
                       vector<double>& best_params, double& best_value) {
    int n = model.params.size();
    vector<double> current_params = best_params;
    double current_value = best_value;
    
    for (int i = 0; i < n; i++) {
        // Пробуем шаг в положительном направлении
        vector<double> params_positive = current_params;
        params_positive[i] += step_sizes[i];
        double f_positive = model.J(params_positive);
        
        if (f_positive < current_value) {
            current_params = params_positive;
            current_value = f_positive;
            continue;  // Переходим к следующей координате
        }
        
        // Пробуем шаг в отрицательном направлении
        vector<double> params_negative = current_params;
        params_negative[i] -= step_sizes[i];
        double f_negative = model.J(params_negative);
        
        if (f_negative < current_value) {
            current_params = params_negative;
            current_value = f_negative;
        }
        // Если оба шага неудачны - остаемся на месте
    }
    
    // Обновляем лучшие найденные параметры
    best_params = current_params;
    best_value = current_value;
}

// Основной метод оптимизации
vector<double> optimize(Model& model, double eps = 1e-10, double init_step = 0.1,
                        double accel = 2.0, double reduction = 0.5, int max_iter = 10000, double delta = 1e-3) {
    int n = model.params.size();
    
    // Текущая точка и значение функции
    vector<double> x_current = model.params;
    double f_current = model.J(x_current);
    
    // Лучшая точка (результат исследующего поиска)
    vector<double> x_base = x_current;
    double f_base = f_current;
    
    // Шаги по каждому направлению
    vector<double> step_sizes(n, init_step);
    for(int i = 0; i < n; i++) {
        step_sizes[i] = model.params[i] * init_step;
    }
    
    int iteration = 0;
    
    while (iteration < max_iter) {
        // Шаг 1: Исследующий поиск вокруг текущей базовой точки
        vector<double> x_new = x_base;
        double f_new = f_base;
        exploratory_search(model, step_sizes, x_new, f_new);
        
        // Проверяем, был ли исследующий поиск успешным
        bool search_successful = (f_new < f_base - eps);
        
        if (search_successful) {
            // Шаг 2: Успешный исследующий поиск - делаем поиск по образцу
            vector<double> x_pattern(n);
            for (int i = 0; i < n; i++) {
                x_pattern[i] = x_new[i] + accel * (x_new[i] - x_base[i]);
            }
            
            double f_pattern = model.J(x_pattern);
            
            // Если поиск по образцу успешен, принимаем новую точку
            if (f_pattern < f_new) {
                x_base = x_new;     // Новая базовая точка
                f_base = f_new;
                x_current = x_pattern;  // Продолжаем из точки образца
                f_current = f_pattern;
            } else {
                // Поиск по образцу неудачен, остаемся в точке исследующего поиска
                x_base = x_new;
                f_base = f_new;
                x_current = x_new;
                f_current = f_new;
            }
        } else {
            // Шаг 3: Исследующий поиск неудачен - уменьшаем шаги
            bool all_steps_small = true;
            for (int i = 0; i < n; i++) {
                if (step_sizes[i] > eps) {
                    step_sizes[i] *= reduction;
                    all_steps_small = false;
                }
            }
            
            // Если все шаги стали очень маленькими - завершаем
            if (all_steps_small && f_base < delta) {
                cout << "Оптимизация завершена на итерации " << iteration << endl;
                break;
            }
            
            // Сбрасываем базовую точку для нового исследующего поиска
            x_current = x_base;
            f_current = f_base;
        }
        
        // Обновляем параметры модели для следующей итерации
        model.params = x_current;
        
        // Вывод информации о текущей итерации (опционально)
        if (iteration % 100 == 0) {
            cout << "Итерация " << iteration << ": f = " << f_current << endl;
            cout << "Параметры: Ron=" << model.params[0] << ", Roff=" << model.params[1] 
                    << ", mu_v=" << model.params[2] << ", Vp=" << model.params[3] 
                    << ", Vn=" << model.params[4] << ", D=" << model.params[5] 
                    << ", x0=" << model.params[6] << endl;
        }
        
        iteration++;
    }
    
    // Возвращаем лучшие найденные параметры
    model.params = x_base;
    return x_base;
}


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

int main() {
    ///////////////////////////////////////////////////////////////// Чтение данных ////////////////////////////////////////////////////////////////////////////

    auto data = readCSV("Pt-HfO2-TiN.csv");
    
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

    int steps = 1000;
    Model model(t_exp, V_exp, I_exp, steps);

    optimize(model);
    
    // Решаем систему уравнений
    cout << "Решаем систему уравнений мемристора..." << endl;

    cout << model.J() << endl;
    cout << model.J(model.params) << endl;
    
    cout << "Решение успешно найдено!" << endl;

    // Сохраняем результаты
    ofstream out_file("memristor_results.csv");
    out_file << "t,V,I,x\n";
    
    for (size_t i = 0; i <= steps; i++) {
        out_file << model.t_grid[i] << "," 
                << model.V_grid[i] << ","
                << model.I_sol[i] << "," 
                << model.x_sol[i] << "\n";
    }
    
    out_file.close();
    
    cout << "Результаты сохранены в memristor_results.csv" << endl;
    cout << "Количество точек расчета: " << model.t_grid.size() << endl;
    

    return 0;
}