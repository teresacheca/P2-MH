#include <iostream>
#include <fstream>
#include <vector>
#include "../include/random.h"
#include <stdlib.h>
#include <algorithm>
#include "../include/BL.h"



#ifndef MEMETICO_H
#define MEMETICO_H

using namespace std;

struct cromosoma_m
    {
        vector<int> genes;      //Vector de genes, tiene el tamaño de la población y está compuesto de 0 y 1
        double valor;           //Es el valor de la función objetivo, servirá para comparar los cromosomas. 
        bool necesita_evaluarse = true;     //Variable que indica si es necesario reevaluar el individuo o no. En un incio se encuentra a true
    };

class Memetico{
    public:

        //n=nº de alelos totales; m=nº de alelos con valor 1
        int n, m;

        //Numero de cromosomas que tendrá la muestra
        int tam = 50;

        //Vector de padres
        std::vector<cromosoma_m> padres;

        //Vector de hijos
        std::vector<cromosoma_m> hijos;

        int evaluaciones = 100000;

        cromosoma_m mejor;

        BL bl;

        //El tipo será 0, 1 o 2, dependiendo de a cuantos cromosomas le aplicamos la búsqueda local
        //  -0: se aplicará a toda la población
        //  -1: se aplicará a un 10% de la población
        //  -2: se aplicará a los cromosomas que se encuentren entre los 10% mejores
        int tipo;


        //Tipo cruce: representa un índice
        //  - si vale 0: El tipo de cruce que utilizaremos será CRUCE UNIFORME
        //  - si vale 1: El tipo de cruce que utilizaremos será CRUCE POR POSICIÓN
        int tipo_cruce;

        std::vector<std::vector<float>> matriz_valores;

        //Constructor, se le pasa el nombre del fichero que contiene los datos para calcular los valores de los cromosomas 
        //y el tamaño de las muestras
        Memetico(string f, int t, int t_c);     

        void leer_matriz(string f);

        void generar_cromosomas();

        double evalua_individuo(vector<int> v);

        double funcion_objetivo(cromosoma_m ind);

        double funcion_principal();

        vector<cromosoma_m> operador_seleccion(vector<cromosoma_m> p);

        cromosoma_m encuentra_mejor(vector<cromosoma_m> v);
       
        cromosoma_m operador_cruce1(cromosoma_m padre1, cromosoma_m padre2);

        vector<int> reparador(vector<int> genes);

        double obtener_contribucion(int e, vector<int> v);

        double genera_contribuciones(vector<int> gen, int i);

        pair<cromosoma_m, cromosoma_m> operador_cruce2(cromosoma_m padre1, cromosoma_m padre2);

        cromosoma_m mutacion(cromosoma_m c);

        vector<cromosoma_m> aplicar_BL(vector<cromosoma_m> c, vector<int> i_aplicar);

        vector<int> poblacion_a_aplicar(vector<cromosoma_m> c);

        vector<int> seleccionar_indice_mejores(double porcentaje);



};

#endif