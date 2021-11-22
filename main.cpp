#include <iostream>
#include "iomanip"

using namespace std;

const int N = 10;

//функция по варианту
double f(double x) {
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

  for (int i = 1; i <= N; i++) {
    y3[i] = 2 * (1 - L2) * y2[i] + L2 * (y2[i + 1] - y2[i - 1]) - y1[i];
  }
  return y3;
}

int main() {

  auto time = clock();

  double t = 0.05;
  double a = 1;
  double b = 0.5;

  double h = a / N;
  double tau = b / N;

  //u0
  auto *y1 = new double[N + 1];
  //u1
  auto *y2 = new double[N + 1];

  auto **u = new double *[N + 1];
  for (int i = 0; i <= N; i++) {
    u[i] = new double[N + 1];
  }

  for (int i = 0; i <= N; i++) {
    y1[i] = funcExample(i * h);
    y2[i] = y1[i] + tau * g(i);
    u[i][0] = y1[i];
    u[i][1] = y2[i];
  }

  auto *y3 = pp(tau, h, y1, y2);
  t += tau;

  for (int k = 2; k <= N; k++) {
    for (int i = 0; i <= N; i++) {
      u[i][k] = y3[i];
      y1[i] = y2[i];
      y2[i] = y3[i];
    }
    y3 = pp(tau, h, y1, y2);
    t = t + tau;
  }

  for (int i = 0; i <= N; i++) {
    for (int j = 0; j <= N; j++) {
      cout << u[i][j] << " ";
    }
    cout << "\n";
  }

  delete[] y1;
  delete[] y2;
  delete[] y3;
  for (int i = 0; i < N; i++) {
    delete[] u[i];
  }
  return 0;
}
