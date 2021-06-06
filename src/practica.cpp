#include <iostream>
#include "../include/generacional.h"
#include "../include/BL.h"
#include <fstream>
#include <vector>
//#include "../include/random.h"
#include <stdlib.h>
#include <algorithm>
#include "../include/timer.h"
#include "../include/estacionario.h"
#include "../include/memetico.h"

float calcula_desviacion(float mejor, float obtenido){
  float desviacion;

  desviacion = abs(100*(mejor-obtenido)/mejor);

  return desviacion;

}


int main(int argc, char *argv[]){

  if (argc <= 1) {
		cerr<<"\nEl programa necesita un fichero y un float como parametro\n\n";
		return 0;
	}

  double valor;
  double time;
  float desviacion;

  float mejor_coste = (float)strtod(argv[2],NULL);

  start_timers();

  //Generacional con cruce uniforme
  Generacional generacional1(argv[1], 0);
  valor = generacional1.funcion_principal();
  time= elapsed_time();
  desviacion = calcula_desviacion(mejor_coste, valor);
  cout << "Tiempo: " << time<< " seg"<< endl;
  cout << "Desviacion: " << desviacion << endl; 

  //Generacinal con cruce basado en posicion
  start_timers();
  Generacional generacional2(argv[1], 1);
  valor = generacional2.funcion_principal();
  time= elapsed_time();
  desviacion = calcula_desviacion(mejor_coste, valor);
  cout << "Tiempo: " << time << " seg"<< endl;
  cout << "Desviacion: " << desviacion << endl; 

  //Estacionario con cruce uniforme
  start_timers();
  Estacionario estacionario1(argv[1], 0);
  valor = estacionario1.funcion_principal();
  time= elapsed_time();
  desviacion = calcula_desviacion(mejor_coste, valor);
  cout << "Tiempo: " << time << " seg"<< endl;
  cout << "Desviacion: " << desviacion << endl; 

  //Estacionario con cruce basado en posición
  start_timers();
  Estacionario estacionario2(argv[1], 1);
  valor = estacionario2.funcion_principal();
  time= elapsed_time();
  desviacion = calcula_desviacion(mejor_coste, valor);
  cout << "Tiempo: " << time << " seg"<< endl;
  cout << "Desviacion: " << desviacion << endl; 

 
  //AM-(10,1.0)
  start_timers();
  Memetico memetico0_1(argv[1], 0, 1);
  valor = memetico0_1.funcion_principal();
  time= elapsed_time();
  desviacion = calcula_desviacion(mejor_coste, valor);
  cout << "Tiempo: " << time << " seg"<< endl;
  cout << "Desviacion: " << desviacion << endl; 


  //AM-(10,0.1)
  start_timers(); 
  Memetico memetico1_1(argv[1], 1, 1);
  valor = memetico1_1.funcion_principal();
  time= elapsed_time();
  desviacion = calcula_desviacion(mejor_coste, valor);
  cout << "Tiempo: " << time << " seg"<< endl;
  cout << "Desviacion: " << desviacion << endl;

  //AM-(10,0.1mej)
  start_timers(); 
  Memetico memetico2_1(argv[1], 2, 1);
  valor =memetico2_1.funcion_principal();
  time= elapsed_time();
  desviacion = calcula_desviacion(mejor_coste, valor);
  cout << "Tiempo: " << time << " seg"<< endl;
  cout << "Desviacion: " << desviacion << endl; 

}