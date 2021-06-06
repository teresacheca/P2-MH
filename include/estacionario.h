#include <iostream>
#include <fstream>
#include <vector>
#include "../include/random.h"
#include <stdlib.h>
#include <algorithm>



#ifndef ESTACIONARIO_H
#define ESTACIONARIO_H

using namespace std;

struct cromosoma_e
    {
        vector<int> genes;      //Vector de genes/alelos, tiene el tamaño de la población y está compuesto de 0 y 1
        double valor;           //Es el valor de la función objetivo, servirá para comparar los cromosomas. 
        bool necesita_evaluarse = true;     //Variable que indica si es necesario reevaluar el individuo o no. En un incio se encuentra a true
    };

class Estacionario{
    public:

        //n=nº de alelos totales; m=nº de alelos con valor 1
        int n, m;

        //Numero de cromosomas que tendrá la muestra
        int tam = 50;

        //Vector de padres
        std::vector<cromosoma_e> padres;

        //Pareja de hijos
        pair<cromosoma_e,cromosoma_e> hijos;

        int evaluaciones = 100000;

       //cromosoma_e mejor;

        
        std::vector<std::vector<float>> matriz_valores;

        //Tipo cruce: representa un índice
        //  - si vale 0: El tipo de cruce que utilizaremos será CRUCE UNIFORME
        //  - si vale 1: El tipo de cruce que utilizaremos será CRUCE POR POSICIÓN
        int tipo_cruce;

        //Constructor, se le pasa el nombre del fichero que contiene los datos para calcular los valores de los cromosomas 
        //y el tamaño de las muestras
        Estacionario(string f, int t_c);     

        void leer_matriz(string f);

        void generar_cromosomas();

        double evalua_individuo(vector<int> v);

        double funcion_objetivo(cromosoma_e ind);

        double funcion_principal();

        pair<cromosoma_e, cromosoma_e> operador_seleccion(vector<cromosoma_e> p);

        cromosoma_e encuentra_mejor(vector<cromosoma_e> v);
       
        cromosoma_e operador_cruce1(cromosoma_e padre1, cromosoma_e padre2);

        vector<int> reparador(vector<int> genes);

        double obtener_contribucion(int e, vector<int> v);

        double genera_contribuciones(vector<int> gen, int i);

        pair<cromosoma_e, cromosoma_e> operador_cruce2(cromosoma_e padre1, cromosoma_e padre2);

        cromosoma_e mutacion(cromosoma_e c);

        pair<int,int> encuentra_peores(vector<cromosoma_e> v);

        cromosoma_e maximo(cromosoma_e c1, cromosoma_e c2, cromosoma_e c3);




};

#endif