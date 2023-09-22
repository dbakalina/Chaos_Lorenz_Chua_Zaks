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

void Lorenz_model(double* point, double* dPoint, double* parameters) {
    //Обозначения:
    //point={x,y,z}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={sigma(0), r(1), beta(2)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    dPoint[0] = -parameters[0] * (point[0] - point[1]);
    dPoint[1] = point[0] * (parameters[1] - point[2]) - point[1];
    dPoint[2] = point[0] * point[1] - parameters[2] * point[2];

}

void PWS_Lorenz_model(double* point, double* dPoint, double* parameters) {
    //Обозначения:
    //point={x,y,z}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={Omega(0),b(1),alpha(2),delta(3),lambda(4),nu(5),omega(6)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    if ((abs(point[0]) < 1) && (point[2] < parameters[1]))
    {
        //система Аs
        dPoint[0] = point[0];
        dPoint[1] = -parameters[2] * point[1];
        dPoint[2] = -parameters[5] * point[2];
    }
    else
    {
        if (((point[0] < -1) && (point[2] <= parameters[1]))
            || ((point[0] < -sign(point[1])) && (point[2] > parameters[1])))
        {
            //система Al
            dPoint[0] = -parameters[4] * (point[0] + 1) + parameters[6] * (point[2] - parameters[1]);
            dPoint[1] = -parameters[3] * (point[1] + parameters[0]);
            dPoint[2] = -parameters[6] * (point[0] + 1) - parameters[4] * (point[2] - parameters[1]);
        }
        else
        {
            //система Ar
            dPoint[0] = -parameters[4] * (point[0] - 1) - parameters[6] * (point[2] - parameters[1]);
            dPoint[1] = -parameters[3] * (point[1] - parameters[0]);
            dPoint[2] = parameters[6] * (point[0] - 1) - parameters[4] * (point[2] - parameters[1]);
        }
    }
}

int main()
{
    setlocale(LC_ALL, "Russian");
    int dim = 3;

    //point={x(0),y(1),z(2)}
    double* point = (double*)malloc(sizeof(double) * dim);
    double* point_pwl = (double*)malloc(sizeof(double) * dim);
    //dPoint={dx/dt, dy/dt, dz/dt}
    double* dPoint = (double*)malloc(sizeof(double) * dim);
    double* dPoint_pwl = (double*)malloc(sizeof(double) * dim);
    //parameters={alpha(0), beta(1), a(2), b(3)}
    double* parameters = (double*)malloc(sizeof(double) * 3);
    double* parameters_pwl = (double*)malloc(sizeof(double) * 7);

    ofstream f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
    f1.open("equilibrium_lorenz.txt", ios::out);
    f2.open("equilibrium_lorenz_pwl.txt", ios::out);

    //parameters={sigma(0), r(1), b(2)}
    parameters[0] = 10.0;     //sigma
    parameters[1] = 24.4;     //r
    parameters[2] = 8.0 / 3.0;  //b

    //parameters={Omega(0),b(1),alpha(2),delta(3),lambda(4),nu(5),omega(6)}
    parameters_pwl[0] = 1;      //Omega
    parameters_pwl[1] = 3.8;    //b........3.8
    parameters_pwl[2] = 2.0;    //alpha
    parameters_pwl[3] = 0.588;  //delta
    parameters_pwl[4] = 0.294;  //lambda
    parameters_pwl[5] = 0.65;   //nu
    parameters_pwl[6] = 2.0;    //omega

    double dt = 0.001;  //Шаг для метода
    int time = 20000;  //Время работы метода

    //начальные значения
    double lambda = -0.5 * (parameters[0] + 1) + sqrt(0.25 * (parameters[0] + 1) * (parameters[0] + 1) + parameters[0] * (parameters[1] - 1));

    point[0] = 0.001;
    point[1] = ((parameters[0] + lambda) / parameters[0]) * 0.001;
    point[2] = 0;

    point_pwl[0] = 0.01;
    point_pwl[1] = 0;
    point_pwl[2] = 0;

    //Состояния равновесия

    const double x = sqrt(parameters[2] * (parameters[1] - 1));
    const double y = sqrt(parameters[2] * (parameters[1] - 1));
    const double z = parameters[1] - 1;

    /*
    const double x_pwl = 1;
    const double y_pwl = 1;
    const double z_pwl = parameters_pwl[1];
    */

    f1 << "0  " << "0  " << "0   " << endl << x << "   " << y << "   " << z << "  " << endl << -x << "  " << -y << "  " << z << endl;
    //f2 << "0  " << "0  " << "0   " << endl << x_pwl << "   " << y_pwl << "   " << z_pwl << "  " << endl << -x_pwl << "  " << -y_pwl << "  " << z_pwl << endl;

    int num = 0;
    cout << "1 - Lorenz model, 2 - PWL_Lorenz" << endl;
    cin >> num;

    if (num == 1) {
        cout << "Lorenz model" << endl;

       f3.open("point_lorenz_24_4.txt", ios::out);
        time = 40000;

        //parameters[0] = 9.67;
        parameters[1] = 24.4;

        //for (int i = 0; i < 3; i++) {
        //    cout << parameters[i] << endl;
        //}

        cout << "sigma = " << parameters[0] << endl;   //sigma
        cout << "r = " << parameters[1] << endl;       //r
        cout << "b = " << parameters[2] << endl;   //b

        for (int i = 0; i < time; i++) {
            RK(Lorenz_model, point, dPoint, parameters, &dim, &dt);
            f3 << i << " " << point[0] << " " << point[1] << " " << point[2] <<
               " " << -point[0] << " " << -point[1] << " " << point[2] << endl;

        }
    }
    if (num == 2) {
        cout << "PWL_Lorenz model" << endl;

        f4.open("point_lorenz_pwl_2.81.txt", ios::out);
        time = 30000;

        //parameters_pwl[5] = 1.1;

        parameters_pwl[5] = 1.25;

        parameters_pwl[1] = 2.81;

        cout << "Omega = " << parameters_pwl[0] << endl;   //Omega
        cout << "b = " << parameters_pwl[1] << endl;       //b
        cout << "alpha = " << parameters_pwl[2] << endl;   //alpha
        cout << "delta = " << parameters_pwl[3] << endl;   //delta
        cout << "lambda = " << parameters_pwl[4] << endl;  //lambda
        cout << "nu = " << parameters_pwl[5] << endl;      //nu
        cout << "omega = " << parameters_pwl[6] << endl;  //omega

        for (int i = 0; i < time; i++) {
            RK(PWS_Lorenz_model, point_pwl, dPoint_pwl, parameters_pwl, &dim, &dt);
            //if (i > 16000) {
            f4 << i << " " << point_pwl[0] << " " << point_pwl[1] << " " << point_pwl[2] <<
                " " << -point_pwl[0] << " " << -point_pwl[1] << " " << point_pwl[2] << endl;
            // }

        }

        const double x_pwl = 1;
        const double y_pwl = 1;
        const double z_pwl = parameters_pwl[1];

        f2 << "0  " << "0  " << "0   " << endl << x_pwl << "   " << y_pwl << "   " << z_pwl << "  " << endl << -x_pwl << "  " << -y_pwl << "  " << z_pwl << endl;
    }

    if (num == 3) {
        cout << "Map x(gamma) for pwl_Lorenz model" << endl;

        f5.open("map_x_gamma_pwl_nu_1.25.txt", ios::out);
        time = 100000;

        parameters_pwl[5] = 1.25;
        //parameters_pwl[5] = 0.65;
        // for (int i = 0; i < 7; i++) {
        //     cout << parameters_pwl[i] << endl;
        // }
        cout << "Omega = " << parameters_pwl[0] << endl;   //Omega
        cout << "b = " << parameters_pwl[1] << endl;       //b
        cout << "alpha = " << parameters_pwl[2] << endl;   //alpha
        cout << "delta = " << parameters_pwl[3] << endl;   //delta
        cout << "lambda = " << parameters_pwl[4] << endl;  //lambda
        cout << "nu = " << parameters_pwl[5] << endl;      //nu
        cout << "omega = " << parameters_pwl[6] << endl;  //omega

        point_pwl[0] = 0.1; //Начальный у равен 0

        double gamma_cr = (exp((-parameters_pwl[4] / parameters_pwl[6]) * atan(parameters_pwl[4] / parameters_pwl[6]))) *
            2 * sqrt(1 + (parameters_pwl[4] * parameters_pwl[4]) / (parameters_pwl[6] * parameters_pwl[6]));

        double r = exp((-3 * 3.14 * parameters_pwl[3]) / (2 * parameters_pwl[6]));


        double gamma = 0;

        double eps = 0.001;

        double tmp_x = point_pwl[0];

        for (gamma = 0.6; gamma <= 2; gamma += 0.001)
        {
            //cout << gamma << endl;
            point_pwl[0] = 0.1;
            //tmp_x = 0.1; //Начальный у равен 0
            // (double j = 0.1; j < 1; j += 0.01) {

                //point_pwl[0] = j;
                //cout << point_pwl[0] << endl;

            for (int i = 0; i < 1000; i++) {

                tmp_x = point_pwl[0];

                if (tmp_x > 0) {
                    //Отображение х
                    point_pwl[0] = 1.0 - gamma + gamma * pow(tmp_x, parameters_pwl[5]);
                }

                else {
                    //Отображение х
                    point_pwl[0] = gamma - 1.0 - gamma * pow(abs(tmp_x), parameters_pwl[5]);
                }

                if (i > 980) {
                    //if (abs(point_pwl[0] - tmp_x) < eps) {
                    f5 << gamma << " " << point_pwl[0] << " " << -point_pwl[0] << endl;
                    //}


                }
                // cout << point_pwl[0] << endl;

            }

        }
    }


    if (num == 4) {
        cout << "Map x(gamma) for pwl_Lorenz model" << endl;

        f6.open("map_x_gamma_pwl_nu_0.65.txt", ios::out);
        
        parameters_pwl[5] = 0.65;

        cout << "Omega = " << parameters_pwl[0] << endl;   //Omega
        cout << "b = " << parameters_pwl[1] << endl;       //b
        cout << "alpha = " << parameters_pwl[2] << endl;   //alpha
        cout << "delta = " << parameters_pwl[3] << endl;   //delta
        cout << "lambda = " << parameters_pwl[4] << endl;  //lambda
        cout << "nu = " << parameters_pwl[5] << endl;      //nu
        cout << "omega = " << parameters_pwl[6] << endl;  //omega

        point_pwl[0] = 0.1; //Начальный x

        double gamma = 0;
        double tmp_x = point_pwl[0];

        for (gamma = 0.5; gamma <= 2; gamma += 0.001)
        {
            cout << gamma << endl;
            point_pwl[0] = 0.1;

            for (int i = 0; i < 900; i++) {

                tmp_x = point_pwl[0];

                if (tmp_x > 0) {
                    //Отображение х
                    point_pwl[0] = 1.0 - gamma + gamma * pow(tmp_x, parameters_pwl[5]);
                }

                else {
                    //Отображение х
                    point_pwl[0] = gamma - 1.0 - gamma * pow(abs(tmp_x), parameters_pwl[5]);
                }

                if (i > 800) {
                    f6 << gamma << " " << point_pwl[0] << " " << -point_pwl[0] << endl;
                }
                // cout << point_pwl[0] << endl;

            }

        }
    }

    if (num == 5) {
        //Построение отображений в координатах 

        f7.open("map_x_x+1_gamma_1.325.txt", ios::out);

        parameters_pwl[5] = 1.25;

        cout << "Omega = " << parameters_pwl[0] << endl;   //Omega
        cout << "b = " << parameters_pwl[1] << endl;       //b
        cout << "alpha = " << parameters_pwl[2] << endl;   //alpha
        cout << "delta = " << parameters_pwl[3] << endl;   //delta
        cout << "lambda = " << parameters_pwl[4] << endl;  //lambda
        cout << "nu = " << parameters_pwl[5] << endl;      //nu
        cout << "omega = " << parameters_pwl[6] << endl;  //omega

        point_pwl[0] = 0.1; //Начальный x

        double gamma = 1.325;
        double tmp_x = point_pwl[0];

        cout << gamma << endl;
        point_pwl[0] = 0.1;

        for (int i = 0; i < 2; i++) {

            tmp_x = point_pwl[0];

            if (tmp_x > 0) {
                //Отображение х
                point_pwl[0] = 1.0 - gamma + gamma * pow(tmp_x, parameters_pwl[5]);
            }

            else {
                //Отображение х
                point_pwl[0] = gamma - 1.0 - gamma * pow(abs(tmp_x), parameters_pwl[5]);
            }

            if (i > 0) {
                f7 << tmp_x << " " << point_pwl[0] << " "<< -tmp_x << " " << -point_pwl[0] << endl;
            }
            // cout << point_pwl[0] << endl;

        }
            


    }

}
