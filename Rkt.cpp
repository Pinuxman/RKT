#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>


using namespace std;


#define PI 3.14159265



int main() {
	std::ofstream out("Protocol.txt", std::ios::app);
	out.setf(ios::left);


	cout.setf(ios::left);
	

	double mr = 298.4 * 1000 + 5700; //масса ракеты
	double C_f = 0.4;//коэффициент лобового сопротивления
	double S = 7;//площадь сечения ракеты

	double ro_0 = 1.225;
	double ro = 1.225; //плотность воздуха 
	double g = 9.8;
	double R_z = 6.371*pow(10,6);	//без 10^6
	double M_z = 5.972*pow(10,24);	//без 10^24
	double G = 6.67*pow(10,-11);		//без 10^-11


	double a_x = 0;
	double v_x = 0;
	double a_y = 0;
	double v_y = 0;
	double h = 0;
	double x = 0;
	double y = R_z;
	double r = R_z;

	double t = 0;
	
	// 1-й этап
	double P_1 = 3506600; //тяга 1
	double P_2 = 887650;  //тяга 2
	double P_3 = 294000;  //тяга 3
	double Ugol_1 = PI/2 - 1.4;
	double Ugol_2 = PI/2 - 0.483;
	double Ugol_3 = PI/2 - 0.0648;
	double fi;

	double dm; //Дельта М 

	double L_1 = 2520; //скорость 1
	double L_2 = 2430; //скорость 2 
	double L_3 = 3200; //скорость 3 
	
	dm = P_1 / L_1 + P_2 / L_2;
	double Fsopr_x = 0;
	double Fsopr_y = 0;
	double dt = 0.01;// дискретный шаг времени

	setlocale(LC_ALL, "rus");
	bool flag_vx = true;

	//первый этап

	for (; t <= 119; t+=dt) {
		ro = ro_0 * (1 - h / 118000); // плотность воздуха
		if (h < 118000) {
			Fsopr_x = C_f * (ro * (pow(v_x, 2) / 2)) * S;
			Fsopr_y = C_f * (ro * (pow(v_y, 2) / 2)) * S;
		}
		else {
			Fsopr_x = 0;
			Fsopr_y = 0;
		}
		
		r = sqrt(pow(x,2)+pow(y,2));
		fi = acos(y/r);
		if (x < 0) {
			fi = 2*PI - fi;
		}
		g = (G * M_z) / pow(r, 2);
		a_x = ((P_1 + P_2) * sin(Ugol_1 + fi) - mr*g*sin(fi) - Fsopr_x) / mr;
		a_y = ((P_1 + P_2) * cos(Ugol_1 + fi) - mr*g*cos(fi) - Fsopr_y) / mr;

		v_x += a_x * dt; // умножаем на дельту t
		mr -= dm * dt; // Масса ракеты
		v_y += a_y * dt;
		x += v_x*dt;
		y += v_y*dt;
		h = r - R_z;
		if (h > 118000) {
			if ( flag_vx){
				v_x += 329*cos(fi);
				v_y += 329*sin(fi);
				flag_vx = false;
			}
		}

	}
	cout  << "t = "  << round(t);
	out << "t = " << round(t);
	cout << "\tmr = " << mr;
	out << "\tmr = " << mr;
	cout << "\tV_x = " << v_x;
	out << "\tV_x = " << v_x;
	cout  << "\ta_x = " << a_x;
	out << "\ta_x = " << a_x;
	cout << "\tV_y = " << v_y;
	out << "\tV_y = " << v_y;
	cout << "\ta_y = " << a_y;
	out << "\ta_y = " << a_y;
	cout << "\th = " << h;
	out << "\th = " << h;
	cout << "\tfi = " << fi*180/PI  << endl;
	out << "\tfi = " << fi*180/PI << endl;
	
	cout << endl << "-----------------------------------------------------------------" << endl << "Второй этап" << endl;
	out << endl << "-----------------------------------------------------------------" << endl << "Второй этап" << endl;
	//второй этап
	
	dm = P_2 / L_2;

	for (t = 120; t <= 302; t+=dt)
	{
		ro = ro_0 * (1 - h / 118000); // плотность воздуха
		if (h < 118000) {
			Fsopr_x = C_f * (ro * (pow(v_x, 2) / 2)) * S;
			Fsopr_y = C_f * (ro * (pow(v_y, 2) / 2)) * S;
		}
		else {
			Fsopr_x = 0;
			Fsopr_y = 0;
		}

		r = sqrt(pow(x,2)+pow(y,2));
		fi = acos(y/r);
		if (x < 0) {
			fi = 2*PI - fi;
		}
		g = (G * M_z) / pow(r, 2);
		a_x = ((P_2) * sin(Ugol_2 + fi) - mr*g*sin(fi) - Fsopr_x) / mr;
		a_y = ((P_2) * cos(Ugol_2 + fi) - mr*g*cos(fi) - Fsopr_y) / mr;

		v_x += a_x * dt; // умножаем на дельту t
		mr -= dm * dt; // Масса ракеты
		v_y += a_y * dt;
		x += v_x*dt;
		y += v_y*dt;
		h = r - R_z;
		if (h > 118000) {
			if (flag_vx) {
				v_x += 329*cos(fi);
				v_y += 329*sin(fi);
				flag_vx = false;
			}
		}
	}

	cout  << "t = "  << round(t);
	out << "t = " << round(t);
	cout << "\tmr = " << mr;
	out << "\tmr = " << mr;
	cout << "\tV_x = " << v_x;
	out << "\tV_x = " << v_x;
	cout  << "\ta_x = " << a_x;
	out << "\ta_x = " << a_x;
	cout << "\tV_y = " << v_y;
	out << "\tV_y = " << v_y;
	cout << "\ta_y = " << a_y;
	out << "\ta_y = " << a_y;
	cout << "\th = " << h;
	out << "\th = " << h;
	cout << "\tfi = " << fi*180/PI  << endl;
	out << "\tfi = " << fi*180/PI << endl;


	cout << endl << "-----------------------------------------------------------------" << endl << "Третий этап" << endl;
	out << endl << "-----------------------------------------------------------------" << endl << "Третий этап" << endl;
	//третий этап
	//return 3;
	dm = P_3 / L_3;

	for (t = 303; t < 543; t+=dt)
	{
		ro = ro_0 * (1 - h / 11800);
		if (h < 118000) {
			Fsopr_x = C_f * (ro * (pow(v_x, 2) / 2)) * S;
			Fsopr_y = C_f * (ro * (pow(v_y, 2) / 2)) * S;
		}
		else {
			Fsopr_x = 0;
			Fsopr_y = 0;
		}
		


		r = sqrt(pow(x,2)+pow(y,2));
		fi = acos(y/r);
		if (x < 0) {
			fi = 2*PI - fi;
		}
		g = (G * M_z) / pow(r, 2);
		a_x = ((P_3) * sin(Ugol_3 + fi) - mr*g*sin(fi) - Fsopr_x) / mr;
		a_y = ((P_3) * cos(Ugol_3 + fi) - mr*g*cos(fi) - Fsopr_y) / mr;

		v_x += a_x * dt; // умножаем на дельту t
		mr -= dm * dt; // Масса ракеты
		v_y += a_y * dt;
		x += v_x*dt;
		y += v_y*dt;
		h = r - R_z;
	}
		cout  << "t = "  << round(t);
		out << "t = " << round(t);
		cout << "\tmr = " << mr;
		out << "\tmr = " << mr;
		cout << "\tV_x = " << v_x;
		out << "\tV_x = " << v_x;
		cout  << "\ta_x = " << a_x;
		out << "\ta_x = " << a_x;
		cout << "\tV_y = " << v_y;
		out << "\tV_y = " << v_y;
		cout << "\ta_y = " << a_y;
		out << "\ta_y = " << a_y;
		cout << "\th = " << h;
		out << "\th = " << h;
		cout << "\tfi = " << fi*180/PI  << endl;
		out << "\tfi = " << fi*180/PI << endl;

		double a;
		double T = 0;
		double orb_v = 0;
		orb_v = sqrt(pow(v_x,2)+pow(v_y,2));
		a = 40*pow(10,13)*(R_z + h)/(80*pow(10,13)-(pow(orb_v,2))*(R_z + h));
		T = (2*PI*sqrt(pow(a,3)/(40*pow(10,13))))/60;

		t = 0;
		int i = 0;
		double apo = 0;
		double per = INFINITY;
		for(;t < T*60;t+=dt){
			i++;
			g = (G * M_z) / pow(r, 2);
			a_x = -g*sin(fi);
			a_y = -g*cos(fi);
			v_x += a_x*dt;
			v_y += a_y*dt;
			x += v_x*dt;
			y += v_y*dt;
			r = sqrt(x*x+y*y);
			h = r - R_z;
			fi = acos(y/r);
			if (x < 0) {
				fi = 2*PI - fi;
			}
			if(i % 5000 == 0){
				cout << " t = " << t;
				cout << "\t h = " << h;
				cout << "\t r = " << r;
				cout << "\t fi = " << fi*180/PI << endl; 
			}
			if(h < per){
				per = h;
			}
			if(h > apo){
				apo = h;
			}
		}
		cout << "apo = " << apo << endl;
		cout << "per = " << per << endl;
		cout << "orb_v = " << orb_v << endl;
		cout << "a = " << a << endl;
		cout << "T = " << T << endl;
		out.close();
	return 0;
}

