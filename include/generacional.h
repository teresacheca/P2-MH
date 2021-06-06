#include <iostream>
#include <fstream>
#include <vector>
#include "../include/random.h"
#include <stdlib.h>
#include <algorithm>



#ifndef GENERACIONAL_H
#define GENERACIONAL_H

using namespace std;

struct cromosoma
    {
        vector<int> genes;      //Vector de genes/alelos, tiene el tamaño de la población y está compuesto de 0 y 1
        double valor;           //Es el valor de la función objetivo, servirá para comparar los cromosoma. 
        bool necesita_evaluarse = true;     //Variable que indica si es necesario reevaluar el individuo o no. En un incio se encuentra a true
    };

class Generacional{
    public:

        //n=nº de alelos totales; m=nº de alelos con valor 1
        int n, m;

        //Numero de cromosoma que tendrá la muestra
        int tam = 50;

        //Vector de padres
        std::vector<cromosoma> padres;

        //Vector de hijos
        std::vector<cromosoma> hijos;

        int evaluaciones = 100000;

        cromosoma mejor;

        
        std::vector<std::vector<float>> matriz_valores;

        //Tipo cruce: representa un índice
        //  - si vale 0: El tipo de cruce que utilizaremos será CRUCE UNIFORME
        //  - si vale 1: El tipo de cruce que utilizaremos será CRUCE POR POSICIÓN
        int tipo_cruce;

        //Constructor, se le pasa el nombre del fichero que contiene los datos para calcular los valores de los cromosoma 
        //y el tamaño de las muestras
        Generacional(string f, int t_c);     

        void leer_matriz(string f);

        void generar_cromosomas();

        double evalua_individuo(vector<int> v);

        double funcion_objetivo(cromosoma ind);

        double funcion_principal();

        vector<cromosoma> operador_seleccion(vector<cromosoma> p);

        cromosoma encuentra_mejor(vector<cromosoma> v);
       
        cromosoma operador_cruce1(cromosoma padre1, cromosoma padre2);

        vector<int> reparador(vector<int> genes);

        double obtener_contribucion(int e, vector<int> v);

        double genera_contribuciones(vector<int> gen, int i);

        pair<cromosoma, cromosoma> operador_cruce2(cromosoma padre1, cromosoma padre2);

        cromosoma mutacion(cromosoma c);



};

#endif