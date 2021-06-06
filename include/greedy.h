
#include <iostream>
#include <vector>
#include "../include/random.h"

#include <fstream>
#include <stdlib.h>
#include <algorithm>

#ifndef GREEDY_H
#define GREEDY_H



using namespace std;

class Greedy{
    public:

        //Elementos seleccionados
        std::vector<int> sel;

        //Elementos seleccionados ordenados
        std::vector<int> ordenados;

        //Matriz que contiene las distancias de los elementos elegidos
        std::vector<std::vector<float>> distancias;

        //Vector que contine la probabilidad de que elijan como el mejor (suma de distancias)
        std::vector<float> probabilidad;

        unsigned long semilla;
        
        //n=nº elementos totales; m=nº muestra
        int n, m;

        string fichero;
    
        Greedy();

        ~Greedy();

        Greedy(string f);

        float funcion_objetivo();

        float calcular_total();
        
        //Devuelve el índice del elemento que tiene la distancia menor
        int calcula_peor(std::vector<float> v);

        //Función que lee el fichero de entrada, inicializa n y m, y crea la matriz con las distancas (solo para los elementos elegidos)
        void inicializa_matriz(string f);

        //Devuelve el índice del elemento que tiene la distancia mayor
        int calcula_mejor( std::vector<float> v);

        void crear_seleccionados();

        void calcula_probabilidad();

        std::vector<int> get_ordenados();

        std::vector<float> get_probabilidad();

        std::vector<std::vector<float>> get_matriz_distancias();

        int get_N();

        int get_M();

        float leer_matriz(int peque, int grande);

        

};

#endif