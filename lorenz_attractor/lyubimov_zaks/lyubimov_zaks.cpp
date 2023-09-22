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

void Lubimov_Zaks_model(double* point, double* dPoint, double* parameters) {
    //Обозначения:
    //point={x,y,z}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={sigma(0), r(1), beta(2), D(3)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    //double D = 0;

    dPoint[0] = -parameters[0] * (point[0] - point[1]) + parameters[0] * point[1] * (point[2] - parameters[1]) * parameters[3];
    dPoint[1] = point[0] * (parameters[1] - point[2]) - point[1];
    dPoint[2] = point[0] * point[1] - parameters[2] * point[2];

}

int main()
{
    setlocale(LC_ALL, "Russian");
    int dim = 3;

    //point={x(0),y(1),z(2)}
    double* point = (double*)malloc(sizeof(double) * dim);
    
    //dPoint={dx/dt, dy/dt, dz/dt}
    double* dPoint = (double*)malloc(sizeof(double) * dim);
    
    //parameters={alpha(0), beta(1), a(2), b(3)}
    double* parameters = (double*)malloc(sizeof(double) * 4);
    
    //parameters={sigma(0), r(1), beta(2), D(3)}
    parameters[0] = 10.0;  //sigma
    parameters[1] = 13.34; ///15.66;  //r......28.0....homoclinic = 13.9270
    parameters[2] = 8.0 / 3.0;  //b
    parameters[3] = 0.054;
    //parameters[3] = 0.054;      //D


    // параметры для измененной системы
    parameters[3] = 0.0584;
    parameters[1] = 15.089;

    ofstream f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;

    cout << "sigma = " << parameters[0] << endl;   //sigma
    cout << "r = " << parameters[1] << endl;       //r
    cout << "b = " << parameters[2] << endl;       //b
    cout << "D = " << parameters[3] << endl;       //D

    double dt = 0.001;  //Шаг для метода
    int time = 1000000;

    int num = 1;

    //double lambda = -0.5 * (parameters[0] + 1) + 0.5 * sqrt(pow((parameters[0] - 1), 2) +
                                                            //4 * parameters[1] * parameters[0] * (1 - parameters[1]));

    point[0] = 0.01;
    point[1] = 0.1;//point[1] = ((parameters[0] + lambda) / parameters[0]) * 0.001;
    point[2] = 0;

    double x = -(sqrt(sqrt(-parameters[2] * parameters[2] * (4* parameters[3] - 1)* parameters[1] * parameters[1])
    + parameters[2]* parameters[1] - 2* parameters[2])/sqrt(2));
    
    double y = (sqrt(sqrt(-parameters[2] * parameters[2] * (4 * parameters[3] - 1) * parameters[1] * parameters[1])
        + parameters[2] * parameters[1] - 2 * parameters[2]) * 
        (sqrt(-parameters[2] * parameters[2] * (4 * parameters[3] - 1) * parameters[1] * parameters[1]) - parameters[2] * parameters[1])) / (2*sqrt(2) * parameters[1] * parameters[2] * parameters[3]);
    
    double z = ((sqrt(-parameters[2] * parameters[2] * (4 * parameters[3] - 1) * parameters[1] * parameters[1])
        + parameters[2] * parameters[1] - 2 * parameters[2]) *
        (-sqrt(-parameters[2] * parameters[2] * (4 * parameters[3] - 1) * parameters[1] * parameters[1]) + parameters[2] * parameters[1])) / (4 * parameters[1] * parameters[2] * parameters[2] * parameters[3]);

    cout << "1 - Point Lyubimov_Zaks, 2 - Map" << endl;
    cin >> num;

    if (num == 1) {
        time = 35000;
        //time = 1E+6;
        f2.open("equilibrium_Lyubimov_Zaks.txt", ios::out);
        f2 << "0  " << "0  " << "0   " << endl << x << "   " << y << "   " << z << "  " << endl << -x << "  " << -y << "  " << z << endl;

      f1.open("lyubimov_H_2_mod.txt", ios::out);
        for (int i = 0; i < time; i++) {
                RK(Lubimov_Zaks_model, point, dPoint, parameters, &dim, &dt);
                //if (i > 900000) {
                f1 << i << " " << point[0] << " " << point[1] << " " << point[2] <<
                    " " << -point[0] << " " << -point[1] << " " << point[2] << endl;
            //}
        }
    }


    if (num == 3) {
        f6.open("equilbrium.txt", ios::out);
        for (parameters[1] = 11; parameters[1] <= 17.5; parameters[1] += 1) {
            x = -(sqrt(sqrt(-parameters[2] * parameters[2] * (4 * parameters[3] - 1) * parameters[1] * parameters[1])
                + parameters[2] * parameters[1] - 2 * parameters[2]) / sqrt(2));
            z = ((sqrt(-parameters[2] * parameters[2] * (4 * parameters[3] - 1) * parameters[1] * parameters[1])
                + parameters[2] * parameters[1] - 2 * parameters[2]) *
                (-sqrt(-parameters[2] * parameters[2] * (4 * parameters[3] - 1) * parameters[1] * parameters[1]) + parameters[2] * parameters[1])) / (4 * parameters[1] * parameters[2] * parameters[2] * parameters[3]);
            f6 << parameters[1] << " " << x << " "<< -x << " " << z << endl;
        }
    }

    if (num == 4) {
        double tmp = point[2];
        int counter = 0;

        double z_rl = 0;
        double x_rl = 0;

        f3.open("equlibrium_x.txt", ios::out);
        f4.open("Map_x_r_new.txt", ios::out);
        f5.open("tmp.txt", ios::out);

        cout << "sigma = " << parameters[0] << endl;   //sigma
        cout << "b = " << parameters[2] << endl;       //b
        cout << "D = " << parameters[3] << endl;       //D

        for (parameters[1] = 11.9; parameters[1] <= 17.5; parameters[1] += 0.01) {
            cout << "r = " << parameters[1] << endl;
            point[0] = 0.01;
            point[1] = 0.1;
            point[2] = 0;

            z_rl = (2 * parameters[3] * parameters[1] + sqrt(1 - 4 * parameters[3]) - 1) / (2 * parameters[3]);

            for (int i = 0; i < time; i++) {
                RK(Lubimov_Zaks_model, point, dPoint, parameters, &dim, &dt);
                
                if (i > 250000) {
                    dt = 0.001;
                    //cout << tmp << endl;
                    //f5 << parameters[1] << " " << point[2] << " " << z_rl << endl;
                    if ((point[2] < z_rl) && (tmp > z_rl)) {
                        x = -(sqrt(sqrt(-parameters[2] * parameters[2] * (4 * parameters[3] - 1) * parameters[1] * parameters[1])
                            + parameters[2] * parameters[1] - 2 * parameters[2]) / sqrt(2));
                        f4 << parameters[1] << " " << point[0] << endl;
                        
                        //cout << tmp << endl;
                        counter++;
                    }
                    //if ((-point[0] > z_rl & -tmp < z_rl)) {
                       // f4 << parameters[1] << " " << -point[0] << endl;
                       // counter++;
                    //}
                    if (counter == 1000) { counter = 0; break; }
                }
                tmp = point[2];
            }
        }
    }


    if (num == 5) {
        double tmp = point[2];
        int counter = 0;

        double z_rl = 0;
;
        f4.open("Map_x_r_new.txt", ios::out);
        //f5.open("tmp.txt", ios::out);

        cout << "sigma = " << parameters[0] << endl;   //sigma
        cout << "b = " << parameters[2] << endl;       //b
        cout << "D = " << parameters[3] << endl;       //D

        for (parameters[1] = 11.5; parameters[1] <= 17.5; parameters[1] += 0.01) {
            cout << "r = " << parameters[1] << endl;
            point[0] = 0.01;
            point[1] = 0.1;
            point[2] = 0;

            z_rl = (2 * parameters[3] * parameters[1] + sqrt(1 - 4 * parameters[3]) - 1) / (2 * parameters[3]);

            for (int i = 0; i < time; i++) {
                RK(Lubimov_Zaks_model, point, dPoint, parameters, &dim, &dt);

                if (i > 250000) {
                    dt = 0.001;
                    if ((point[2] < z_rl) && (tmp > z_rl)) {
                        //f4 << parameters[1] << " " << point[0] << endl;
                        counter++;
                    }
                    
                    if (counter > 200) {
                        f4 << parameters[1] << " " << point[0] << endl;
                    }
                    if (counter == 1000) {
                        counter = 0;
                        break;
                    }
                }
                tmp = point[2];
            }
        }
    }

    if (num == 6) {
        time = 1E+6;
        double tmp = point[2];
        int counter = 0;

        double z_rl = 0;

        f4.open("Map_x_r_modify.txt", ios::out);

        cout << "sigma = " << parameters[0] << endl;   //sigma
        cout << "b = " << parameters[2] << endl;       //b
        //cout << "D = " << parameters[3] << endl;       //D

        //Меняем параметр D
        for (parameters[3] = 0.095; parameters[3] > 0.04; parameters[3] -= 0.0001) {
            cout << "D = " << parameters[3] << endl;

            parameters[1] = 0.881 / parameters[3];  //r = 0.0881 / D
            cout << "r = " << parameters[1] << endl;

            point[0] = 0.1;
            point[1] = 0.1;
            point[2] = 0;

            z_rl = (2 * parameters[3] * parameters[1] + sqrt(1 - 4 * parameters[3]) - 1) / (2 * parameters[3]);

            for (int i = 0; i < time; i++) {
                RK(Lubimov_Zaks_model, point, dPoint, parameters, &dim, &dt);

                if (i > 700000) {
                    //cout << point[2] << endl;
                    if ((point[2] < z_rl) && (tmp > z_rl)) {
                        f4 << parameters[3] << " "<< parameters[1] << " " << point[0] << endl;
                        //cout << parameters[1] << " " << point[0] << endl;
                        //cout << "Work" << endl;
                        counter++;
                    }
                    if (counter == 1000) {
                        counter = 0;
                        break;
                    }
                }
                tmp = point[2];
            }
        }
    }

    if (num == 7) {
        for (parameters[3] = 0.093; parameters[3] > 0.04; parameters[3] -= 0.0001) {
            

            parameters[1] = 0.881 / parameters[3];  //r = 0.0881 / D
            
        }

    }
}

