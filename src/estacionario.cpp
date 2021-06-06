#include <iostream>
#include "../include/estacionario.h"
#include <fstream>


Estacionario::Estacionario(string f, int t_c){

    tipo_cruce = t_c;

    //Abrimos el fichero
    ifstream ifile;
    ifile.open(f);

    if(ifile.is_open() == 0){
        cout << "NO ESTA ABIERTO" << endl;
    }

    //Leemos la primero línea que contiene el valor de n y m, número total de alelos y número de alelos con valor 1 respectivamente
    double linea;
    ifile >> linea;
    n = linea;
    ifile >> linea;
    m=linea;

    cout << "N: " << n << "   M:  " << m << endl;

    leer_matriz(f);

    //Generamos el grupo de cromosomas inicial que vamos a utilizar
    generar_cromosomas();

    
}

//Leemos la matriz de datos que nos servirá para saber el valor de cada individuo/cromosoma
void Estacionario::leer_matriz(string f){

    ifstream ifile(f);
    std::vector<float> aux;
    float fila, columna, valor;
    ifile >> fila;
    ifile >> columna;

    matriz_valores.resize(n);

    for( int i=0 ; i<n ; i++){
        matriz_valores[i].resize(n);
    }
  
   
    while(!ifile.eof()){
        
        ifile >> fila;
        ifile >> columna;
        ifile >> valor;
        
        matriz_valores[fila][columna] = valor;
        matriz_valores[columna][fila] = valor;
        
    }
}

//Generar_cromosomas genera el primer frupo de cromosomas de forma aleatoria dentro de la población total de la que disponemos 
void Estacionario::generar_cromosomas(){

    //Crearemos un individuo
    cromosoma_e ind; 

    //En el vector de genes de esta individuo añadiremos primero el número de 0's necesarios y luego el número de 1's necesarios
    for(int k =0; k<m; k++){
        ind.genes.push_back(1);
    }
    for(int k =0; k<n-m; k++){
        ind.genes.push_back(0);
    }

    //Como queremos tam individuos, que tengan el mismo número de 0's y 1's que tiene el inidviduo creado
    //Haremos un bucle de tam iteraciones que mezclará el vector de genes del inidviduo, evaluará al individuo y lo añadirá al vecotr de padres
    for(int i=0; i<tam; i++){
        std::random_shuffle(ind.genes.begin(), ind.genes.end());
        ind.valor = funcion_objetivo(ind);
        ind.necesita_evaluarse = false;
        padres.push_back(ind);
    }
    
}


//Esta función recorre el vectro de cromosomas del individuo. Luego, guarda en un vector los índices de los alelos que tienen valor 1.
//Cuando tenemos el vecotr de índices, llamamos a la función objetivo, que nos devolverá el valor correspondicente el vetor de índices
double Estacionario::funcion_objetivo(cromosoma_e ind){
    vector<int> indices;
    vector<int>::iterator it;
    it = ind.genes.begin();
    int contador =0;
    for(it; it!=ind.genes.end(); it++){
        if((*it) == 1){
            indices.push_back(contador);
        }

        contador++;
    }

    return evalua_individuo(indices);
}

//Función objetivo: recorre el valor de ínidices pasado como parámetro,
//hace la sumatoria de los valores de los ínidices con el resto
double Estacionario::evalua_individuo(vector<int> v){
    double suma=0;
    vector<int>::iterator it1, it2;
    it1 = it2 = v.begin();
    for(it1; it1 != v.end()-1; it1++){
        it2 = it1+1;
        for( it2; it2 != v.end(); it2++){
            suma += matriz_valores[(*it1)][(*it2)];
        }
    }

    return suma;
}

//Es la función que realiza todas las funciones necesarias para el algoritmo
double Estacionario::funcion_principal(){
    
    int contador =0;
    while(contador < evaluaciones){
        //Ya tenemos los padres generados

        //SELECCIÓN
        //Seleccionamos a los dos individuos que actuarán de padres con el operador de selección     
        hijos = operador_seleccion(padres);
  
        //CRUCE
        //Hacemos el cruce entre los dos padres aplicando los operadores de cruce
        //En este caso, la probabilidad de cruca será siempre 1, por lo que siempre generaremos los dos hijos
        switch (tipo_cruce)
        {
            case 0:
                //Primer operador de cruce: Cruce uniforme
                //Como queremos obtener dos hijos y el operador sólo genera uno, lanzamos el operados dos veces con los mismos padres
                hijos.first = operador_cruce1(hijos.first, hijos.second);
                hijos.second = operador_cruce1(hijos.first, hijos.second);
            break;
            
            case 1:
                //Segundo operador de cruce: Cruce basado en posición
                //Genera una pareja de hijos
                pair<cromosoma_e, cromosoma_e> h = operador_cruce2(hijos.first, hijos.second);
                hijos.first = h.first;
                hijos.second= h.second;
                    
            break;
        }
        

        //MUTACIÓN
        //Aplicamos el operador de mutación
           
        //Generamos un número aleatorio entre 0 y 1. 
        //Si es menor de 0.1/num. alelos, el cromosoma_e mutará, si no, el cromosoma_e se mantendrá igual
        double probabilidad_mutacion = ((double) rand() / (RAND_MAX));
        if(probabilidad_mutacion <= (0.1/n)){
            hijos.first = mutacion(hijos.first);
        }
        probabilidad_mutacion = ((double) rand() / (RAND_MAX));
        if(probabilidad_mutacion <= (0.1/n)){
            hijos.second = mutacion(hijos.second);
        }


        //EVALUACIÓN Y REEMPLAZO
        //Los hijos obtenidos compiten por entrar en el grupo de padres
        //Evaluamos los hijos de nuevo y los comparamos con los peores de la población,
        //el que tenga mayor valor, entrará en el vector de padres
        hijos.first.valor = funcion_objetivo(hijos.first);
        hijos.second.valor = funcion_objetivo(hijos.second);

       
        pair<int,int> peores = encuentra_peores(padres);
        pair<cromosoma_e, cromosoma_e> padres_peores;
        padres_peores.first = padres.at(peores.first);
        padres_peores.second = padres.at(peores.second);


        cromosoma_e auxiliar;
        cromosoma_e y1 = maximo(hijos.first, hijos.second, auxiliar);
        cromosoma_e y2 = maximo(padres_peores.first, padres_peores.second, auxiliar);
        cromosoma_e mejor1 = maximo(y1, y2, auxiliar);

        y1 = maximo(hijos.first, hijos.second, mejor1);
        y2 = maximo(padres_peores.first, padres_peores.second, mejor1);
        cromosoma_e mejor2 = maximo(y1, y2, mejor1);

        padres.at(peores.first) = mejor1;
        padres.at(peores.second) = mejor2;
        
        vector<cromosoma_e>::iterator it;
        it = padres.begin();
        for(it; it!=padres.end(); it++){
            contador++;
            if((*it).necesita_evaluarse == true){
                (*it).valor = funcion_objetivo((*it));
                (*it).necesita_evaluarse = false;

            }
        }

    }
    cromosoma_e mejor = encuentra_mejor(padres);
    return mejor.valor;

}

//Operador de selección: de un grupo de cromosomas, eligirá dos aleatorios, estos serán los padres seleccionados
pair<cromosoma_e, cromosoma_e> Estacionario::operador_seleccion(vector<cromosoma_e> p){
    pair<cromosoma_e, cromosoma_e> h;
    int indice1, indice2;
    cromosoma_e individuo1, individuo2;
    

    //Generamos dos números aleatorios entre 0 y 49. 
    //Estos corresponderán a los índices de los cromosomas seleccionados
    indice1 = rand()%50;
    indice2 = rand()%50;

    //Guardamos los cromosomas que se encuentran en los ínidices 
    //generados anteriormente de forma aleatoria
    individuo1 = p.at(indice1);
    individuo2 = p.at(indice2);

    //Ahora añadiremos ambos cromosomas al vector que devolverá la función
    h.first = individuo1;
    h.second = individuo2;


    //Se devuelve la pareja de padres generada
    return h;

}


//Encuentra mejor: Recorre el vector de cromosomas y devuelve el que tenga mayor valor
cromosoma_e Estacionario::encuentra_mejor(vector<cromosoma_e> v){
    vector<cromosoma_e>::iterator it;
    it = v.begin();
    cromosoma_e m;
    m= (*it);
    for(it; it!=v.end(); it++){
        if((*it).valor > m.valor){
            m = (*it);
        }
    }

    return m;
}


//Operador de cruce 1 (Operador uniforme): De dos padres dados, se quiere generar un hijo
//Para ello, primero, recorre el vector de genes de los padres, 
//si el valor del gen coincide, se añade al vector de genes del hijo a generar
//si no coincide, se añade un 0 o 1 inidistintamente
//Una vez hecho esto, aplica el reparador al vector de genes del hijo
cromosoma_e Estacionario::operador_cruce1(cromosoma_e padre1, cromosoma_e padre2){

    vector<int> genes_hijos;
    vector<int>::iterator it1, it2;
    it1 = padre1.genes.begin();
    it2 = padre2.genes.begin();

    //Se recorre el vector de genes de los padres y se comprueban que los valores son iguales
    for(it1; it1 != padre1.genes.end(); it1++){
        //Si los valores son iguales, se guarda este valor en el vector de genes del hijo
        if((*it1) == (*it2)){
            genes_hijos.push_back((*it1));
        }else{
            //Si es distinto, se añade un 0 o un 1 indistintamente
            genes_hijos.push_back((rand() % 2)); //0 o 1 generado que se guarda en un arreglo
        }
    }

    //Se aplica el reparador a los genes evalua_individuol hijo
    genes_hijos = reparador(genes_hijos);
    cromosoma_e hijo;
    hijo.genes= genes_hijos;

    //Necesita evaluarse del hijo se pone a true, para saber que se ha modificado el cromosoma
    //hijo.necesita_evaluarse = true;


    return hijo;

}


//Reparador: Recorre un vector de genes y contabiliza el número de 1's que tiene:
//   -si tuviese el número de 1's inidicado, no lo modifica
//   -si tiene más 1's de los necesarios cambia el valor de los 1's que se encuentran en las posiciones que 
//   generan mayor contribución al valor del inivididuo hasta tener el número de 1's indicado
//   -si tiene menos 1's de los necesariios, cambiará el valor de los 0's que se encuentren en las posiciones
//   que generen mayor contribución al valor del individuo hasta tener el número de 1's indicado
vector<int> Estacionario::reparador(vector<int> genes){
 
    //Primero se recorre el vector de genes y se contabilizan los 1's
    vector<int>::iterator it;
    it = genes.begin();
    int contador_unos = 0;
    for(it; it!=genes.end(); it++){
        if((*it) == 1){
            contador_unos++;
        }
    }

    //Si el número de 1's vale m, devolvemos el vector igual
    if(contador_unos == m){
        return genes;
    }else if(contador_unos > m){
        //Si el número de unos el mayor que m, haremos las iteraciones necesarias hasta que contador_unos sea igual a m
        while (contador_unos > m){
            //Primero, obtendremos las contribuciones de cada alelo del vector total de genes
            //Como solo nos interesan las posiciones en las que tenemos 1's, 
            //las contribuciones en las que el valor del gen sea 0, estas valdrán 0 (aunque no sea su valor real)
            it = genes.begin();
            vector<double> contribuciones;
            int contador =0;
            contribuciones.clear();
            for(it; it!=genes.end(); it++){
                if((*it) == 1){
                    contribuciones.push_back(genera_contribuciones(genes, contador));
                }else{
                    contribuciones.push_back(0);
                }
                contador++;
            }

            //Comprobaremos las contribuciones, y guardaremos la mayor
            vector<double>::iterator it2;
            it2 = contribuciones.begin();
            double mayor=0;
            int indice =0;
            contador =0;
            for(it2; it2 != contribuciones.end(); it2++){
                if((*it2) > mayor){
                    mayor = (*it2);
                    indice = contador;
                }
                contador ++;
            }

            //Cambiaremos el 1 que se encuentra en la posición con mayor contribución por un 0
            genes.at(indice) = 0;

            //restaremos un 1 al contador de unos, puesto que hemos cambiado el valor de unos de ellos
            contador_unos--;
        }


    }else if(contador_unos < m){
        //Si el número de unos es menor a m, haremos las iteraciones necesarias hasta que contador_unos sea igual a m
        while (contador_unos < m){

            //evalua_individuos contribuciones de cada gen del vector de genes
            //Como solo nos interesan las posiciones en las que tenemos 0's, 
            //las contribuciones en las que el valor del gen sea 1, estas valdrán 0 (aunque no sea su valor real)
            it = genes.begin();
            vector<double> contribuciones;
            int contador =0;
            contribuciones.clear();
            for(it; it!=genes.end(); it++){
                if((*it) == 0){
                    contribuciones.push_back(genera_contribuciones(genes, contador));
                }else{
                    contribuciones.push_back(0);
                }
                contador++;
            }

           //Comprobamos las contribuciones y guardamos que que tenga mayor valor
            vector<double>::iterator it2;
            it2 = contribuciones.begin();
            double mayor=0;
            int indice =0;
            contador =0;
            for(it2; it2 != contribuciones.end(); it2++){
                if((*it2) > mayor){
                    mayor = (*it2);
                    indice = contador;
                }
                contador ++;
            }

           //Cambiamos el 0 que se encuentra en la posición de mayor contribución por un 1
            genes.at(indice) = 1;

            //Añdimos un uno al contador de unos
            contador_unos++;


        }
    }

    return genes;
}


//Obtener contribución: devuelve la sumatoria de los valores del punto indicado con el resto de los puntos del vector
double Estacionario::obtener_contribucion(int e, vector<int> v){
    double contribucion = 0;
    vector<int>::iterator it;
    it = v.begin();

    for(it; it != v.end(); it++){
        if((*it) != e){
            contribucion += matriz_valores[(*it)][e];
        }
    }
    return contribucion;
}


//Genera_contribuciones: como sólo nos interesan los alelos que tienen valor 1 del array pasado como parámetro
//la función recorrerá el array y guardará los índices de los alelos que tienen valor 1
//Luego llamará a la funcion obtener_contribucion para el ínidice pasado como paŕametro y el vecotr de índices generado
double Estacionario::genera_contribuciones(vector<int> gen, int i){

    vector<int> indices;
    vector<int>::iterator it;
    it = gen.begin();
    int contador =0;
    for(it; it!=gen.end(); it++){
        if((*it) == 1){
            indices.push_back(contador);
        }

        contador++;
    }

   return  obtener_contribucion(i, indices);

}


//Segundo operador de cruce: OPERADOR BASADO EN POSICIÓN
//El operador obtendrá dos hijos de los dos padres pasados como parámetro
//Para ello, primero se crearán dos hijos llamados "hijo1" e "hijo2"
//Los vectores de genes de dichos hijos se inicializarán a -1
//Luego, se recorrerán los vectores de genes de los padres:
//  - Si los valores los alelos de estos coninciden, este valor se añadirá al vector de genes 
//  de los dos hijos en la posición encontrada
//  - Si los valores no coinciden, se guardaran en un vector de "restos" de los alelos (uno para cada padre)
//Cuando se haya hecho esto, se mezclarán los dos vectores de genes restantes de los padres y:
//  -Al primer hijo, se le irán añadiendo los valores que se encuentren en el vector de restos del padre 1 en las
//posiciones en las que no se han añadido alelos en el hijo (estas tendrán valor -1)
//  -Al segundo hijo, se le irán añadiendo los valores que se encuentren en el vector de restos del padre 2 en las
//posiciones en las que no se han añadido alelos en el hijo (estas tendrán valor -1)
pair<cromosoma_e, cromosoma_e> Estacionario::operador_cruce2(cromosoma_e padre1, cromosoma_e padre2){
    pair<cromosoma_e, cromosoma_e> hijos_generados;
    cromosoma_e hijo1, hijo2;

    //Se inicializan el vector de genes de los hijos
    for(int i=0; i<padre1.genes.size(); i++){
        hijo1.genes.push_back(-1);
        hijo2.genes.push_back(-1);
    }

    vector<int> restos_padre1, restos_padre2;

    //se recorren los vectores de genes de los padres dados
    for(int i=0; i<padre1.genes.size(); i++){
        //Si los valores coinciden, se añaden en la misma posición en los vectores de genes de los hijos
        if(padre1.genes.at(i) == padre2.genes.at(i)){
            hijo1.genes.at(i) = padre1.genes.at(i);
            hijo2.genes.at(i) = padre1.genes.at(i);
        }else{
            //Si no coinciden, se guardan dichos valores en los vectores de restos para cada padre respectivamente
            restos_padre1.push_back(padre1.genes.at(i));
            restos_padre2.push_back(padre2.genes.at(i));
        }
    }

    //Se mezclan los valores de los restos de los genes no seleccionados de los padres
    std::random_shuffle(restos_padre1.begin(), restos_padre1.end());
    std::random_shuffle(restos_padre2.begin(), restos_padre2.end());

    //Añadiremos los genes del padre 1 que no se habían seleccionado al hijo 1
    vector<int>::iterator it;
    it = restos_padre1.begin();
    for(int i=0; i<hijo1.genes.size(); i++){
        if(hijo1.genes.at(i) == -1){
            hijo1.genes.at(i) = (*it);
            it++;
        }
    }

    //Añadiremos los genes del padre 2 que no se habían seleccionado al hijo 2
    it = restos_padre2.begin();
    for(int i=0; i<hijo2.genes.size(); i++){
        if(hijo2.genes.at(i) == -1){
            hijo2.genes.at(i) = (*it);
            it++;
        }
    }

    //Marcamos los indicadores de evaluación de estos hijos a true para aber que es necesario reevaluarlos
  //  hijo1.necesita_evaluarse = true;
   // hijo2.necesita_evaluarse = true;

    hijos_generados.first = hijo1;
    hijos_generados.second = hijo2;

    
    //Devolvemos los hijos generados
    return hijos_generados;

}


//Mutación
//Este algoritmo selecciona un gen aleatorio del vector de genes del cromosoma
//Luego selecciona otro aleatorio:
//  -Si tienen el mismo valor, genera otro de forma aleatoria
//  -Si este tiene distinto valor que el anterior, los intercambia
cromosoma_e Estacionario::mutacion(cromosoma_e c){
    int gen_mutado1, valor1;
    int gen_mutado2, valor2;
    bool salir = false;

    //Generamos un índice de aleatorio entre 0 y 49:
    gen_mutado1 = rand()%50;
    //Guardamos el valor del gen que se encuentra en el índice generado
    valor1 = c.genes.at(gen_mutado1);

    //Generamos índices aleatorios hasta encontrar uno que tenga valor distinto al anteriormente generado
    while(!salir){
        //Generamos un índice de aleatorio entre 0 y 49:
        gen_mutado2 = rand()%50;
        //Guardamos el valor de gen que se encuentra en el ínidice generado
        valor2 = c.genes.at(gen_mutado2);

        //Si los valores son distintos, los intercambiamos y salirmos del bucle
        if(valor1 != valor2){
            c.genes.at(gen_mutado1)= valor2;
            c.genes.at(gen_mutado2) = valor1;
            salir =true;
        }
    }

    //Devolvemos el cromosoma con la mutación del gen
    return c;   

}


pair<int,int> Estacionario::encuentra_peores(vector<cromosoma_e> v){
    vector<cromosoma_e> auxiliar = v;
    vector<cromosoma_e>::iterator it;
    pair<int, int> indices;

    int indice1, indice2;
    double menor;

    
    indice1 =0;
    menor = auxiliar.at(0).valor;
    for(int i=0; i<auxiliar.size(); i++){
        if(auxiliar.at(i).valor < menor){
            indice1 = i;
            menor = auxiliar.at(i).valor;
        }
    }

    indice2 =0;
    menor = auxiliar.at(0).valor;
    for(int i=0; i<auxiliar.size(); i++){
        if(auxiliar.at(i).valor < menor && i!=indice1){
            indice2 = i;
            menor = auxiliar.at(i).valor;
        }
    }

    indices.first = indice1;
    indices.second = indice2;

    return indices;


}

cromosoma_e Estacionario::maximo(cromosoma_e c1, cromosoma_e c2, cromosoma_e c3){
    cromosoma_e mejor;

    if(c1.valor > c2.valor){
        mejor = c1;
        if(c1.valor == c3.valor){
            mejor = c2;
        }
    }else{
        mejor = c2;
         if(c2.valor == c3.valor){
            mejor = c1;
        }
    }

    return mejor;
}

