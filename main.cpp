#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const int N = 10;

//функция по варианту
double func(double x) {
    if (x >= 0.7)
        return (9.0 * x - 4.5) / 0.2;
    if (x < 0.7)
        return (-9.0 * x + 8.1) / 0.2;
    return 0.0;
}

double g(double x) { return 0.0; }

//с методички
double funcExample(double x) {
    if (x < 0.5) return 2.0 * x;
    return 2.0 - 2.0 * x;
}

//Подпрограмма по методичке
double *pp(double tau, double h, const double *y1, const double *y2) {
    //лямбда
    double L = tau / h;
    auto *y3 = new double[N + 1];
    y3[0] = 0;
    y3[N] = 0;
    double L2 = pow(L, 2);

    for (int i = 1; i < N; i++) {
        y3[i] = 2 * (1 - L2) * y2[i] + L2 * (y2[i + 1] - y2[i - 1]) - y1[i];
    }
    return y3;
}

int main() {

    ofstream fileOutput = ofstream("output.txt");
    auto time1 = clock();

    double t = 0.0;
    double a = 1;

    double h = a / N;
    double tau = h * h / 2 / 2;
    double T = 1.0;
    //Находим кол-во точек по шагам в итоге при Т=1 будет 400 точек t
    auto Nt = T / tau;

    //u0
    auto *y1 = new double[N + 1];
    //u1
    auto *y2 = new double[N + 1];

    auto **u = new double *[(long) Nt];

    //начальные и краевые условия
    for (int i = 0; i <= N; i++) {
        auto x = i * h;
        //x от 0 до 1
        y1[i] = func(i * h);
        y2[i] = y1[i] + tau * g(i);
        //u[i][0] = y1[i];
        //u[i][1] = y2[i];
    }

    /*auto *y3 = pp(tau, h, y1, y2);
    t += tau;
    auto row0 = u[0] = new double[N + 1];
    for (int i = 0; i <= N; i++) {
        row0[i] = y3[i];
    }*/
    auto *y3 = new double[N + 1];

    //Зафиксиурем индекс для матрицы вывода
    int uIdx = 0;
    while (t < T) {

        //Перекладывание массивов
        for (int i = 0; i <= N; i++) {
            y1[i] = y2[i];
            y2[i] = y3[i];
        }

        //Вычисление УМФ
        y3 = pp(tau, h, y1, y2);

        int iterIdx = floor(t / T * Nt);
        t = t + tau;

        //перекладывем результат в матрицу вывода. По факту u[i][j] = y3[i]
        if (iterIdx == uIdx || t >= T) {
            double *row = u[uIdx] = new double[N + 1];
            for (int i = 0; i <= N; i++) {
                row[i] = y3[i];
            }
            uIdx++;
        }
    }

    double time2 = (static_cast<double>(clock()) - time1) / CLOCKS_PER_SEC;
    fileOutput << "Time: " << time2 << "\n" << "\n";


    double x = 0.0;
    for (int i = 0; i <= N; i++) {
        fileOutput << x << ",";
        x += 0.1;
    }
    fileOutput << "\n";
    int tCount = 1;
    for (int i = 0; i < uIdx; i++) {
        fileOutput << tCount << ",";
        double *row = u[i];
        for (int j = 0; j <= N; j++) {
            if (j == N) {
                fileOutput << row[j];
            } else
                fileOutput << row[j] << ",";
        }
        tCount++;
        fileOutput << "\n";
    }

    delete[] y1;
    delete[] y2;
    delete[] y3;
    delete[] u;
    return 0;
}
