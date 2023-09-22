#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>

using namespace std;

template<typename T>
int sign(T var)
{
    if (var < 0)
    {
        return -1;
    }
    else { return 1; }
}

//интегратор Рунге-Кутты 4-го порядка с фиксированным шагом
void RK(void(*fun)(double* point, double* dPoint, double* parameters),
    double* point, double* dPoint, double* parameters, int* dim, double* t_step)
{
    //Векторы "точки" длинной с вектор начальных условий
    //vector<double> point(initCon.size());
    //vector<double> buf_point(initCon.size());

    double* buf_point = new double[*dim];

    //Двумерный массив коэффициентов интегрирования
    double** K = new double* [*dim];
    for (int i = 0; i < *dim; i++)
    {
        K[i] = new double[4];
    }

    /////////////////////////////////////////
    //K1
    fun(point, dPoint, parameters);
    for (int i = 0; i < *dim; i++)
    {
        K[i][0] = *t_step * dPoint[i];
        buf_point[i] = point[i] + 0.5 * K[i][0];
    }

    //K2
    fun(buf_point, dPoint, parameters);
    for (int i = 0; i < *dim; i++)
    {
        K[i][1] = *t_step * dPoint[i];
        buf_point[i] = point[i] + 0.5 * K[i][1];
    }

    //K3
    fun(buf_point, dPoint, parameters);
    for (int i = 0; i < *dim; i++)
    {
        K[i][2] = *t_step * dPoint[i];
        buf_point[i] = point[i] + K[i][2];
    }

    //K4
    fun(buf_point, dPoint, parameters);
    for (int i = 0; i < *dim; i++)
    {
        K[i][3] = *t_step * dPoint[i];
    }

    //интегрирование
    for (int i = 0; i < *dim; i++)
    {
        point[i] += (K[i][0] + 2 * K[i][1] + 2 * K[i][2] + K[i][3]) / 6;
    }
    /////////////////////////////////////////

    //очистка памяти
    for (int i = 0; i < *dim; i++)
    {
        delete[] K[i];
    }
    delete[] K;
    delete[] buf_point;

}

void Shimizu_Morioka_model_mod(double* point, double* dPoint, double* parameters) {
    //Обозначения:
    //point={x,y,u}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={alpha(0), lambda(1), gamma(2)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    dPoint[0] = point[1];
    dPoint[1] = point[0] - parameters[1] * point[1] - (point[0] * point[2]) / parameters[2];
    dPoint[2] = -parameters[0] * point[2] + point[0] * point[0];

}

void Shimizu_Morioka_model(double* point, double* dPoint, double* parameters) {
    //Обозначения:
    //point={x,y,z}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={alpha(0), lambda(1)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    dPoint[0] = point[1];
    dPoint[1] = point[0] - parameters[1] * point[1] - point[0] * point[2];
    dPoint[2] = -parameters[0] * point[2] + point[0] * point[0];

}

int main()
{
    setlocale(LC_ALL, "Russian");
    int dim = 3;

    //point={x(0),y(1),z(2)}
    double* point = (double*)malloc(sizeof(double) * dim);

    //dPoint={dx/dt, dy/dt, dz/dt}
    double* dPoint = (double*)malloc(sizeof(double) * dim);

    //parameters={alpha(0), beta(1), mu(2), delta(3)}
    double* parameters = (double*)malloc(sizeof(double) * 2);

    //parameters={alpha(0), beta(1), mu(2), delta(3)}
    parameters[0] = 0.6161;       //alpha
    parameters[1] = 1.0067;     //lambda
    //parameters[2] = 10;     //gamma

    ofstream f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;

    cout << "alpha = " << parameters[0] << endl;      //alpha
    cout << "lambda = " << parameters[1] << endl;     //lambda
    //cout << "gamma = " << parameters[2] << endl;     //gamma

    double dt = 0.001;  //Шаг для метода
    int time = 60000;

    double v = 0.5 * (parameters[1] + sqrt(4 + parameters[1] * parameters[1]));

    point[0] = 0.01;
    point[1] = 0.01;    
    point[2] = 0;

    const double x = sqrt(parameters[0]);
    const double y = 0;
    const double z = 1;

    f1.open("equilibrium_Shimizu_Morioka.txt", ios::out);
    f1 << "0  " << "0  " << "0   " << endl << x << "   " << y << "   " << z << "  " << endl << -x << "  " << y << "  " << z << endl;
    
    f2.open("Shimizu_Morioka_alpha_0.9090_lambda_0.1909.txt", ios::out);
    for (int i = 0; i < time; i++) {
        RK(Shimizu_Morioka_model, point, dPoint, parameters, &dim, &dt);
        if (i > 500) {
            f2 << i << " " << point[0] << " " << point[1] << " " << point[2] <<
                " " << -point[0] << " " << -point[1] << " " << point[2] << endl;
        }

    }
}
