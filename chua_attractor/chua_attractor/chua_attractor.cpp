#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <functional>
#include<fstream>
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

void Chua_model_linear(double* point, double* dPoint, double* parameters)
{
    //Обозначения:
    //point={x(0),y(1),z(2)}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={alpha(0), beta(1), a(2), b(3)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    
    if (point[0] >= 1) {
        //система Аs
        dPoint[0] = parameters[0] * (point[1] - (parameters[3] * point[0] + parameters[2] - parameters[3]));
        dPoint[1] = point[0] - point[1] + point[2];
        dPoint[2] = -parameters[1] * point[1];
    }
    if (abs(point[0]) <= 1) {
        dPoint[0] = parameters[0] * (point[1] - (parameters[2] * point[0]));
        dPoint[1] = point[0] - point[1] + point[2];
        dPoint[2] = -parameters[1] * point[1];
    }
    if (point[0] <= -1) {
        dPoint[0] = parameters[0] * (point[1] - (parameters[3] * point[0] - parameters[2] + parameters[3]));
        dPoint[1] = point[0] - point[1] + point[2];
        dPoint[2] = -parameters[1] * point[1];
    }

    /*

    dPoint[0] = parameters[0] * (point[1] - (0.0625 * pow(point[0], 3) - 0.167 * point[0]));
    dPoint[1] = point[0] - point[1] + point[2];
    dPoint[2] = -parameters[1] * point[1];
    */
}

void Chua_model(double* point, double* dPoint, double* parameters)
{
    //Обозначения:
    //point={x(0),y(1),z(2)}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={alpha(0), beta(1), a(2), c(3)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    dPoint[0] = parameters[0] * (point[1] - (1.0/16.0 * pow(point[0], 3) - 1.0/6.0 * point[0]));
    dPoint[1] = point[0] - point[1] + point[2];
    dPoint[2] = -parameters[1] * point[1];


   // dPoint[0] = parameters[0] * (point[1] + 0.42857 * pow(point[0], 3) + 0.42857 * point[0]);
}

void Return_map(double* point, double* dPoint, double* parameters)
{
    //Обозначения:
    //point={x(0),y(1),z(2)}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={nu(0), omega(1), alpha(2), h(3), Omega(4), beta(5), r(6)}

    double mu = 1 - (1 / parameters[3]) * exp((-3 * 3.14 * parameters[2]) / (2 * parameters[4]));
    double q = exp((-3 * 3.14 * parameters[5]) / (2 * parameters[4]));
    double tmp = 0;

    for (int i = 0; i < 1; i++) {
        tmp = point[0];
        point[0] = mu * sign(point[0]) + (1 - mu) * pow(abs(point[0]), parameters[0]) *
            (point[1] * sin(parameters[1] * log(abs(point[0]))) + point[2] * cos(parameters[1] * log(abs(point[0]))));
        point[1] = q * pow(abs(tmp), parameters[0]) * (point[1] * cos( parameters[1] * log(abs(tmp)) - 
            point[2] * sin(parameters[1] * log(abs(tmp)))));
        point[2] = sign(tmp);
    }
}

void PWL_DS_model(double* point, double* dPoint, double* parameters)
{
    //Обозначения:
    //point={x,y,z}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={nu(0), omega(1), alpha(2), h(3), Omega(4), beta(5), r(6)}

    //вектор производных {dx/dt, dy/dt, dz/dt}

    if ((abs(point[0]) < parameters[3]) && ((point[1] * point[1] + point[2] * point[2] < parameters[6] * parameters[6]) && abs(point[2] < 1)))
    {
        //система Аs
        dPoint[0] = point[0];
        dPoint[1] = -parameters[0] * point[1] + parameters[1] * point[2];
        dPoint[2] = -parameters[1] * point[1] - parameters[0] * point[2];
    }
    else
    {
        if (point[2] > -sign(point[0])) {
            //система Ar
            dPoint[0] = -parameters[2] * (point[0] - parameters[3]) - parameters[4] * (point[2] - 1);
            dPoint[1] = -parameters[5] * point[1];
            dPoint[2] = parameters[4] * (point[0] - parameters[3]) - parameters[2] * (point[2] - 1);
        }
        else {
            //система Al
            dPoint[0] = -parameters[2] * (point[0] + parameters[3]) - parameters[4] * (point[2] + 1);
            dPoint[1] = -parameters[5] * point[1];
            dPoint[2] = parameters[4] * (point[0] + parameters[3]) - parameters[2] * (point[2] + 1);
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
    double* point_pwl_new = (double*)malloc(sizeof(double) * dim);
    double* point_biff = (double*)malloc(sizeof(double) * dim);
    //dPoint={dx/dt, dy/dt, dz/dt}
    double* dPoint = (double*)malloc(sizeof(double) * dim);
    double* dPoint_pwl = (double*)malloc(sizeof(double) * dim);
    //parameters={alpha(0), beta(1), a(2), b(3)}
    double* parameters = (double*)malloc(sizeof(double) * 2);
    double* param = (double*)malloc(sizeof(double) * 7);

    ofstream f, f2,f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
    f.open("point_pwl_0.40512.txt", ios::out);
    f2.open("chua_cubic_9.67.txt", ios::out);
    //f3.open("homoclinic_chua_2.txt", ios::out);
    f4.open("equilibrium_state.txt", ios::out);
    f5.open("map_nu_0. 0.3668.txt", ios::out);
    f6.open("diagram_chua_less_point.txt", ios::out);
    f8.open("diagram_chua__less_point_2.txt", ios::out);
    f7.open("chua_alpha.txt", ios::out);
    f9.open("extended_diagram_PWS_DL.txt", ios::out);
    f10.open("extended_diagram_PWS_DL_2.txt", ios::out);
    f11.open("extended_diagram_PWS_DL_3.txt", ios::out);
    f12.open("equilibrium_pwl.txt", ios::out);


    //parameters={alpha(0), beta(1), a(2), b(3)}
    parameters[0] = 7.1; //10.91865;//6;  //7;(homoclinic)  // alpha  .... 15.6 /// 8.875   //parameters[0] = 8.9;
    parameters[1] = 14.0;  //8.6;(homoclinic)..........//beta    .......25.58  /// 14      homoclinic: 14
   // parameters[2] = -1.0/7.0;//...........a 
  //  parameters[3] = 2.0 / 7.0;    //parameters[3] = 2.0/7.0 - 0.0059992953; - homoclinic......\\b

    //parameters={nu(0), omega(1), alpha(2), h(3), Omega(4), beta(5), r(6)}
    param[0] = 0.75;  //nu.........0.75
    param[1] = 3;     //omega........3
    param[2] = -0.1;  //alpha .......-0.1
    param[3] = 1.60198; // homoclinic
   // param[3] = 1.69; //0.172883;  //h......1.69
    param[4] = 1;  //Omega......1
    param[5] = 20;  //beta..........20
    param[6] = 1.01539;  //r

    double dt = 0.001;  //Шаг для метода
   // int time = 1E6+500000;  //Время работы метода
     int time = 1E6+200000;
   // int time = 2E6;
    //int time = 120000;
   // int time = 100000;

    //начальные значения для chua model

    /*
    //начальные значения для chua
    point[0] = (parameters[3] - parameters[2]) / (parameters[3] + 1) + 0.01; //0.5
    point[1] = 0; //0.1
    point[2] = -(parameters[3] - parameters[2]) / (parameters[3] + 1)+ 0.01;     //
    */

    //начальные значения для chua
    /*
    point[0] = -0.019007; //homoclinic
    point[1] = -0.00220433; //homoclinic
    point[2] = 0.01;
    */

      point[0] = -0.0138288; //0.5
      point[1] = -0.00133472; //0.1
      point[2] = 0.01;


   // point[0] = 0.001; //0.5
   // point[1] = 0.001; //0.1
  //  point[2] = 0.001;

    //начальные значения для pwl_ds
    /*point_pwl[0] = param[3] + 10; //0.5
    point_pwl[1] = 0; //0.1
    point_pwl[2] = 1;
    */

    point_pwl[0] = 0.01; //0.5
    point_pwl[1] = 0; //0.1
    point_pwl[2] = 0;

    //Состояния равновесия для системы Чуа
   // const double x = 1 - parameters[2]/parameters[3];
   // const double y = 0;
   // const double z = parameters[2] / parameters[3] -1;

    const double x = -2*sqrt(2.0/3.0);
    const double y = 0;
    const double z = 2 * sqrt(2.0 / 3.0);

    const double x_pwl = param[3];
    const double y_pwl = 0;
    const double z_pwl = 1;

    f4 << "0  " << "0  " << "0   " << endl << x << "   " << y << "   " << z << "  " << endl << -x << "  " << y << "  " << -z << endl;
    f12 << "0  " << "0  " << "0   " << endl << x_pwl << "   " << y_pwl << "   " << z_pwl << "  " << endl << -x_pwl << "  " << y_pwl << "  " << -z_pwl << endl;

    double tmp = point[0];
    int counter = 0;

    int num = 0;
    cout << "1-Chua model, 2 - PWL_DS, 4 - Diagram, 5 - Exxtended diagram for PWL" << endl;
    cin >> num;


    if (num == 1) {
        cout << "Chua model" << endl;
        point[0] = -0.0148564; //0.5
        point[1] = -0.00144274; //0.1
        point[2] = 0.01;
        //point[0] = -2.7;
        //point[1] = -2;
        //point[2] = 0.01;
        parameters[0] = 9.67;
        parameters[1] = 14;
        //parameters[1] = 11.31484;//гомоклиника
        for (int i = 0; i < 2; i++) {
            cout << parameters[i] << endl;
        }
        time = 39000;
        for (int i = 0; i < time; i++) {
            RK(Chua_model, point, dPoint, parameters, &dim, &dt);
            if (i > 16000) {
                f2 << i << " " << point[0] << " " << point[1] << " " << point[2] <<
                    " " << -point[0] << " " << -point[1] << " " << -point[2] << endl;
            }
            
        }
    }

    if (num == 9) { //эксперимент для кубической Чуа с другими параметрами нелинейности
       // parameters[0] = -4.9261084;
       // parameters[1] = -3.6496350;
        cout << "Diagram" << endl;
        for (int i = 0; i < 2; i++) {
            cout << parameters[i] << endl;
        }
        for (parameters[0] = 8; parameters[0] < 14; parameters[0] += 0.01) {
            cout << "alpha=" << parameters[0] << endl;
            point[0] = -0.0267057; //0.5
            point[1] = -0.00311564; //0.1
            point[2] = 0.01;
            for (int i = 0; i < time; i++) {
                RK(Chua_model, point, dPoint, parameters, &dim, &dt);
                if (i > 300000) {
                    if ((point[1] > 0 & tmp < 0)) {
                        f6 << parameters[0] << " " << point[0] << endl;
                        counter++;
                    }
                    if ((-point[1] > 0 & -tmp < 0)) {
                        f8 << parameters[0] << " " << -point[0] << endl;
                    }
                    //tmp = point[1];
                    if (counter == 1000) { counter = 0; break; }
                }
                tmp = point[1];
            }
        }
       /* cout << "Chua model" << endl;
        for (int i = 0; i < 2; i++) {
            cout << parameters[i] << endl;
        }
        for (int i = 0; i < time; i++) {
            RK(Chua_model, point, dPoint, parameters, &dim, &dt);
            f2 << i << " " << point[2] << " " << point[1] << " " << point[0] <<
                " " << -point[2] << " " << -point[1] << " " << -point[0] << endl;
        }*/
    }

    if (num == 2) {
        param[0] = 0.40512; // nu
        time = 100000;
        cout << "PWL_DS" << endl;
        for (int i = 0; i < time; i++) {
            RK(PWL_DS_model, point_pwl, dPoint_pwl, param, &dim, &dt);
            //if (i > 1900000) {
                f << i << " " << point_pwl[2] << " " << point_pwl[1] << " " << point_pwl[0] <<
                    " " << -point_pwl[2] << " " << -point_pwl[1] << " " << -point_pwl[0] << endl;
            //}
        }
    }

    /*
    if (num == 3) {  //Фигня, которую можно удалить
        cout << "Diagram" << endl;
        for (parameters[3] = (2.0 / 7.0 - 0.005999);  parameters[3] < 1;  parameters[3]+= 0.001) {
            cout << parameters[3] << endl;
            for (double i = 0.000; i < 5.0; i += t) {
                cout << i << " "<< point_biff[0] <<  endl;
                line(point_biff, i); //задание координат секущей для некоторого t
                for (int i = 0; i < time; i++) {
                    RK(Chua_model, point, dPoint, parameters, &dim, &dt);
                    if (point_biff[0] - point[0] < eps && point_biff[1] - point[1] < eps)
                        f6 << parameters[3] << " " << point_biff[0]  << " " << point[0] << " " << point[1] << " " << point[2] << endl;
                    cout << point_biff[0] << " " << point[0] << endl;
                }
            }
        }
    }*/



    if (num == 3) {
       
        //В системе Чуа меняем параметр b
        cout << "Diagram" << endl;
        for (parameters[3] = 2.0 / 7.0 - 0.00599; parameters[3] < 2; parameters[3] += 0.01) {
            cout << parameters[3] << endl;
            for (int i = 0; i < time; i++) {
                RK(Chua_model, point, dPoint, parameters, &dim, &dt);
                if (i > 100000) {
                    if ((point[1] > 0 & tmp < 0)) {
                        f6 << parameters[3] << " " << point[0] << endl;
                        counter++;
                        //cout << counter << " " << (0.1578 * point[0] + 2.01921 * point[1] + 0.15782 * point[2]) << " " << parameters[3] << endl;
                    }
                    //tmp = point[1];
                    if (counter == 1000) { counter = 0; break; }
                }
                tmp = point[1];
            }
        }
    }

    if (num == 4) {
        //В системе Чуа меняем параметр alpha
        cout << "Diagram" << endl;
       // for (parameters[0] = 5; parameters[0] < 20; parameters[0] += 0.01) {
        point[0] = -0.0241628;//-0.0118666; //0.5
        point[1] = -0.00284346;//-0.000851513; //0.1
        point[2] = 0.01;
         for (parameters[0] = 14; parameters[0] > 4; parameters[0] -= 0.01) {
            cout <<"alpha=" << parameters[0] << endl;
            for (int i = 0; i < time; i++) {
                RK(Chua_model, point, dPoint, parameters, &dim, &dt);
                if (i > 300000) {
                    if ((point[1] > 0 & tmp < 0)) {
                        f6 << parameters[0] << " " << point[0] << endl;
                        counter++;
                    }
                    if ((-point[1] > 0 & -tmp < 0)) {
                        f8 << parameters[0] << " " << -point[0] << endl;
                    }
                    //tmp = point[1];
                    if (counter == 1000) { counter = 0; break; }
                }
                tmp = point[1];
            }
        }
    }
    if (num == 14) {
        time = 1E6;
        //В системе Чуа меняем параметр alpha
        cout << "Diagram" << endl;
        // for (parameters[0] = 5; parameters[0] < 20; parameters[0] += 0.01) {
        point[0] = -0.0285906;//-0.0118666; //0.5
        point[1] = -0.00363556;//-0.000851513; //0.1
        point[2] = 0.01;

        parameters[0] = 9; //гомоклиника
        parameters[1] = 11.31484;//гомоклиника
        //parameters[1] = 14;
        cout << "beta=" << parameters[1] << endl;

        for (parameters[0] = 5; parameters[0] < 9.5; parameters[0] += 0.01) {
            cout << "alpha=" << parameters[0] << endl;

            //Начальные значения
            point[0] = -0.0137737;//-0.0118666; //0.5
            point[1] = -0.00143706;//-0.000851513; //0.1
            point[2] = 0.01;

            for (int i = 0; i < time; i++) {
                RK(Chua_model, point, dPoint, parameters, &dim, &dt);
                if (i > 300000) {
                    if ((point[1] > 0 & tmp < 0)) {
                        f6 << parameters[0] << " " << point[0] << endl;
                        //cout << "alpha=" << point[0] << endl;
                        counter++;
                    }
                    if ((-point[1] > 0 & -tmp < 0)) {
                        f8 << parameters[0] << " " << -point[0] << endl;
                    }
                    //tmp = point[1];
                    if (counter == 50) { counter = 0; break; }
                }
                tmp = point[1];
            }
        }
    }

    if (num == 24) {
        //В системе Чуа меняем параметр alpha
        cout << "Diagram" << endl;
        // for (parameters[0] = 5; parameters[0] < 20; parameters[0] += 0.01) {
        point[0] = -0.0285906;//-0.0118666; //0.5
        point[1] = -0.00363556;//-0.000851513; //0.1
        point[2] = 0.01;

        parameters[0] = 9; //гомоклиника
        parameters[1] = 11.31484;//гомоклиника
        cout << "beta=" << parameters[1] << endl;

        for (parameters[0] = 5; parameters[0] < 10; parameters[0] += 0.01) {
            cout << "alpha=" << parameters[0] << endl;
            parameters[1] = 1.3979 * parameters[0] - 1.2185; //обновление значения beta
            //Начальные значения
            point[0] = -0.0137737;//-0.0118666; //0.5
            point[1] = -0.00143706;//-0.000851513; //0.1
            point[2] = 0.01;

            for (int i = 0; i < time; i++) {
                RK(Chua_model, point, dPoint, parameters, &dim, &dt);
                if (i > 300000) {
                    if ((point[1] > 0 & tmp < 0)) {
                        f6 << parameters[0] << " " << point[0] << endl;
                        counter++;
                    }
                    if ((-point[1] > 0 & -tmp < 0)) {
                        f8 << parameters[0] << " " << -point[0] << endl;
                    }
                    //tmp = point[1];
                    if (counter == 700) { counter = 0; break; }
                }
                tmp = point[1];
            }
        }
    }

    if(num == 5) {

        point_pwl[0] = 0.1; //в пределах от 0 до 1
        point_pwl[1] = 0; //Начальный у равен 0
        point_pwl[2] = 1;
        //рассширенная диаграмма для PWS_DL
        //double mu = 1 - (1 / parameters[3]) * exp((-3 * 3.14 * parameters[2]) / (2 * parameters[4]));
        double q = exp((-3 * 3.14 * param[5]) / (2 * param[4]));
        //q = 0;
        double tmp_x = point_pwl[0];
        double tmp_y = point_pwl[1];
        double tmp_z = point_pwl[2];

        double mu = 0;
        param[0] = 0.3668; //Параметр nu
       // param[0] = 0.36718;
       // param[0] = 0.40512;
        param[5] = 160;

        
        //for (mu = 0.1; mu >= -0.15; mu -= 0.0005)
        for (mu = -0.01; mu < 0.01; mu += 0.00001) {
            //cout << mu << endl;
            point_pwl[0] = 0.1; //в пределах от 0 до 1
            point_pwl[1] = 0; //Начальный у равен 0
            point_pwl[2] = 1;

            for (int i = 0; i < 1000; i++) {

                tmp_x = point_pwl[0];
                tmp_y = point_pwl[1];
                tmp_z = point_pwl[2];
                //point[0] = mu * sign(point[0]) + (1 - mu) * pow(abs(point[0]), parameters[0]) * point[2] * cos(parameters[1] * log(abs(point[0])));

                //Отображение х
                point_pwl[0] = mu * sign(tmp_x) + (1 - mu) * pow(abs(tmp_x), param[0]) *
                    (tmp_y * sin(param[1] * log(abs(tmp_x))) + tmp_z * cos(param[1] * log(abs(tmp_x))));
                //Отображение y
                point_pwl[1] = q * pow(abs(tmp_x), param[0]) * (tmp_y * cos(param[1] * log(abs(tmp_x))) - tmp_z * sin(param[1] * log(abs(tmp_x))));
                //Отображение z
                point_pwl[2] = sign(tmp_x);

                if (i > 950) {
                    if (point_pwl[2] <= 1 && point_pwl[2] >= 0) {
                        f5 << mu << " " << point_pwl[0] << endl;
                    }
                }


            }
        }
        
    }

    if (num == 10) {

        point_pwl[0] = 0.1; //в пределах от 0 до 1
        point_pwl[1] = 0; //Начальный у равен 0
        point_pwl[2] = 1;
        //рассширенная диаграмма для PWS_DL
        //double mu = 1 - (1 / parameters[3]) * exp((-3 * 3.14 * parameters[2]) / (2 * parameters[4]));
        double q = exp((-3 * 3.14 * param[5]) / (2 * param[4]));
        q = 0;
        double tmp_x = point_pwl[0];
        double tmp_y = point_pwl[1];
        double tmp_z = point_pwl[2];

        double mu = 0;
        param[0] = 1; //Параметр nu....0.923
        param[5] = 160;

        for (param[0] = 0.25; param[0] < 1; param[0] += 0.001)
            //for (mu = -0.015; mu < 0.15; mu += 0.0001)
        {
            cout << param[0] << endl;
            //point_pwl[0] = 0.1; //в пределах от 0 до 1
           // point_pwl[1] = 0; //Начальный у равен 0
            //point_pwl[2] = 1;

            for (int i = 0; i < 1000; i++) {

                tmp_x = point_pwl[0];
                tmp_y = point_pwl[1];
                tmp_z = point_pwl[2];
                //point[0] = mu * sign(point[0]) + (1 - mu) * pow(abs(point[0]), parameters[0]) * point[2] * cos(parameters[1] * log(abs(point[0])));

                //Отображение х
                point_pwl[0] = mu * sign(tmp_x) + (1 - mu) * pow(abs(tmp_x), param[0]) *
                    (tmp_y * sin(param[1] * log(abs(tmp_x))) + tmp_z * cos(param[1] * log(abs(tmp_x))));
                //Отображение y
                point_pwl[1] = q * pow(abs(tmp_x), param[0]) * (tmp_y * cos(param[1] * log(abs(tmp_x))) - tmp_z * sin(param[1] * log(abs(tmp_x))));
                //Отображение z
                point_pwl[2] = sign(tmp_x);

                if (i > 600) {
                    if (point_pwl[2] <= 1 && point_pwl[2] >= 0) {
                        f5 << param[0] << " " << point_pwl[0] << endl;
                    }
                }


            }
        }

        /*
        for (mu = 0.1; mu >= -0.15; mu -= 0.0005)
        //for (mu = -0.015; mu < 0.15; mu += 0.0001)
        {
            cout << mu << endl;
            //point_pwl[0] = 0.1; //в пределах от 0 до 1
           // point_pwl[1] = 0; //Начальный у равен 0
            //point_pwl[2] = 1;

            for (int i = 0; i < 1000; i++) {

                    tmp_x = point_pwl[0];
                    tmp_y = point_pwl[1];
                    tmp_z = point_pwl[2];
                    //point[0] = mu * sign(point[0]) + (1 - mu) * pow(abs(point[0]), parameters[0]) * point[2] * cos(parameters[1] * log(abs(point[0])));

                    //Отображение х
                    point_pwl[0] = mu * sign(tmp_x) + (1 - mu) * pow(abs(tmp_x), param[0]) *
                        (tmp_y * sin(param[1] * log(abs(tmp_x))) + tmp_z * cos(param[1] * log(abs(tmp_x))));
                    //Отображение y
                    point_pwl[1] = q * pow(abs(tmp_x), param[0]) * (tmp_y * cos(param[1] * log(abs(tmp_x))) - tmp_z * sin(param[1] * log(abs(tmp_x))));
                    //Отображение z
                    point_pwl[2] = sign(tmp_x);

                    if (i > 600) {
                        if (point_pwl[2] <= 1 && point_pwl[2] >= 0) {
                            f5 << mu << " " << point_pwl[0] << endl;
                        }
                    }


            }
        }*/

    }

    if (num == 6) {
        //рассширенная диаграмма для PWS_DL

        double tmp_pwl = point_pwl[0];

        //начальные значения
        //point_pwl[0] = 0.01;
        //point_pwl[1] = 0;
        //point_pwl[2] = 0;

        cout << "parametrs:" << endl;
        for (int i = 0; i < 7; i++) {
            cout << param[i] << endl;
        }
        
        cout << "diagram" << endl;
        for (param[3] = 1.1; param[3] < 6; param[3] += 0.01)
        {
            cout << param[3] << endl;
            point_pwl[0] = 0.01;
            point_pwl[1] = 0;
            point_pwl[2] = 0;
            for (int i = 0; i < time; i++) {
                RK(PWL_DS_model, point_pwl, dPoint_pwl, param, &dim, &dt);
                if (i > 100000) {
                    if (point_pwl[2] < 1 & tmp_pwl > 1) {
                        f10 << param[3] << " " << point_pwl[0] << endl;
                        counter++;
                    }
                    if (counter == 1000) { counter = 0; break; }

                }
                tmp_pwl = point_pwl[2];
            }
        }
        
        /*
        double eps = 0.001;

        cout << "diagram" << endl;
        for (param[3] = 1.1; param[3] < 6; param[3] += 0.01)
        {
            cout << param[3] << endl;

            point_pwl_new[0] = param[3] + eps;
            point_pwl_new[1] = 0;
            point_pwl_new[2] = 1;

            for (int i = 0; i < time; i++) {
                RK(PWL_DS_model, point_pwl_new, dPoint_pwl, param, &dim, &dt);
                if (i > 400000) {
                    if (point_pwl_new[2] < 1 & tmp_pwl > 1) {
                        f11 << param[3] << " " << point_pwl_new[0] << endl;
                        counter++;
                    }
                    if (counter == 1000) { counter = 0; break; }

                }
                tmp_pwl = point_pwl_new[2];
            }
        }
        */
    }

    if (num == 7) {
        double tmp_pwl = point_pwl[0];

        //начальные значения
        //point_pwl[0] = 0.01;
        //point_pwl[1] = 0;
        //point_pwl[2] = 0;

        cout << "parametrs:" << endl;
        for (int i = 0; i < 7; i++) {
            cout << param[i] << endl;
        }

        cout << "diagram" << endl;
        for (param[3] = 1.35; param[3] < 1.4; param[3] += 0.0001)
        {
            cout << param[3] << endl;
            point_pwl[0] = 0.01;
            point_pwl[1] = 0;
            point_pwl[2] = 0;
            for (int i = 0; i < time; i++) {
                RK(PWL_DS_model, point_pwl, dPoint_pwl, param, &dim, &dt);
                if (i > 1200000) {
                    if (point_pwl[2] < 1 & tmp_pwl > 1) {
                        f12 << param[3] << " " << point_pwl[0] << endl;
                        counter++;
                    }
                    if (counter == 1000) { counter = 0; break; }

                }
                tmp_pwl = point_pwl[2];
            }
        }
    }

    if (num == 8) {
        double tmp_pwl = point_pwl[0];

        //начальные значения
        point_pwl[0] = 0.01;
        point_pwl[1] = 0;
        point_pwl[2] = 0;

        cout << "parametrs:" << endl;
        for (int i = 0; i < 7; i++) {
            cout << param[i] << endl;
        }

        cout << "diagram" << endl;
        for (param[3] = 1.4; param[3] > 1.35; param[3] -= 0.0001)
        {
            cout << param[3] << endl;
            //point_pwl[0] = 0.01;
            //point_pwl[1] = 0;
            //point_pwl[2] = 0;
            for (int i = 0; i < time; i++) {
                RK(PWL_DS_model, point_pwl, dPoint_pwl, param, &dim, &dt);
                if (i > 500000) {
                    if (point_pwl[2] < 1 & tmp_pwl > 1) {
                        f12 << param[3] << " " << point_pwl[0] << endl;
                        counter++;
                    }
                    if (counter == 1000) { counter = 0; break; }

                }
                tmp_pwl = point_pwl[2];
            }
        }
    }
}

