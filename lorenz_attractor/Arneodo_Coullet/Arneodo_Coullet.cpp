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

void Arneodo_Coullet_model(double* point, double* dPoint, double* parameters) {
    //Обозначения:
    //point={x,y,z}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={alpha(0), beta(1), mu(2), delta(3)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    //double D = 0;

    dPoint[0] = parameters[0] * (point[0] - point[1]) ;
    dPoint[1] = -4 * parameters[0] * point[1] + point[0] * point[2] + parameters[2] * pow(point[0], 3);
    dPoint[2] = -parameters[3] * parameters[0] * point[2] + point[0] * point[1] + parameters[1] * pow(point[2], 2);

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
    double* parameters = (double*)malloc(sizeof(double) * 4);

    //parameters={alpha(0), beta(1), mu(2), delta(3)}
    parameters[0] = 1.5;       //alpha
    parameters[1] = -0.07;       //beta
    parameters[2] = 0.07;  //mu
    parameters[3] = 0.95;      //delta

    ofstream f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
     
    cout << "alpha = " << parameters[0] << endl;      //alpha
    cout << "beta = " << parameters[1] << endl;       //beta
    cout << "mu = " << parameters[2] << endl;         //mu
    cout << "delta = " << parameters[3] << endl;      //delta

    f3.open("point_Arneodo_Coullet_alpha_2.txt", ios::out);
    double dt = 0.001;  //Шаг для метода
    int time = 50000;


    //double lambda = 0.5 * (sqrt(-4 * parameters[1] * parameters[1] * parameters[0] * parameters[3] + 4 * parameters[1] * parameters[0] + parameters[0] * parameters[0]
      //  - 2 * parameters[0] + 1) - parameters[0] - 1);

    point[0] = 0.01;
    point[1] = 0.1;    //point[1] = ((parameters[0] + lambda) / parameters[0]) * 0.001;
    point[2] = 0;

    const double x = -2 * sqrt(-parameters[0] * parameters[0] * (4 * parameters[1] - parameters[3]));
    const double y = -2 * sqrt(-parameters[0] * parameters[0] * (4 * parameters[1] - parameters[3]));
    const double z = 4 * parameters[0];

    f1.open("equilibrium_Arneodo.txt", ios::out);
    f1 << "0  " << "0  " << "0   " << endl << x << "   " << y << "   " << z << "  " << endl << -x << "  " << -y << "  " << z << endl;
        //<< "0  " << "0  " << "   " << (parameters[0] * parameters[3])/ parameters[1] << endl;

    for (int i = 0; i < time; i++) {
        RK(Arneodo_Coullet_model, point, dPoint, parameters, &dim, &dt);
        if (i > 16000) {
        f3 << i << " " << point[0] << " " << point[1] << " " << point[2] <<
            " " << -point[0] << " " << -point[1] << " " << point[2] << endl;
        }

    }
}


