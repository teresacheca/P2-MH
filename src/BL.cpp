#include <iostream>
#include "../include/BL.h"
#include <fstream>

BL::BL(){
}

void BL::inicializa(std::vector<std::vector<float>> matriz, int _n, int _m){
    n = _n;
    m= _m;
    setMatriz(matriz);

   /* for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

           cout <<  matriz_distancias[i][j] << " ";
        }
        cout << " --" <<endl;
    }*/

    
}

void BL::setMatriz( std::vector<std::vector<float>> matriz){
    matriz_distancias.resize(n);

    for( int i=0 ; i<n ; i++){
        matriz_distancias[i].resize(n);
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            matriz_distancias[i][j] = matriz[i][j];
         //   cout <<  matriz_distancias[i][j] << " ";
        }
       // cout << " --" <<endl;
    }
}

double BL::funcion_principal_bl(){
    int contador =0;
    while(contador < evaluaciones && (explorados.size()< n-m) ){
        vector<puntos>::iterator it;
        it = sel.begin();
        puntos vecino = genera_vecino(sel, explorados);
        puntos menor= *it;
        sel.erase(it);
        
       
        
        double contribucion_vecino = funcion_objetivo(vecino.indice, sel);
        contador ++;
        
        vecino.contribucion = contribucion_vecino;
       
       if(vecino.contribucion > menor.contribucion){
            factorizar(menor.indice, vecino.indice);
            sel.push_back(vecino);
            sel = ordenar_por_contribucion(sel);
            explorados.push_back(menor);
            //generar_contribuciones();
        }else{
            sel.push_back(menor);
            sel = ordenar_por_contribucion(sel);
            explorados.push_back(vecino);
        }

        
    }

    double total = funcion_objetivo_total(sel);
    vector<puntos>::iterator it2;
  /*  cout << "INDICES BL" <<endl;
        it2 = sel.begin();
        for(it2; it2!=sel.end(); it2++){
            cout << (*it2).indice << " ";
        }
        cout << "END " <<endl;*/

    //cout << "Total: " << total <<endl;

    return total;

}

void BL::leer_matriz(string f){
// cout << "Inicializando matriz... "<< endl;
    ifstream ifile(f);
    std::vector<float> aux;
    float fila, columna, distancia;
    ifile >> fila;
    ifile >> columna;

    matriz_distancias.resize(n);

    for( int i=0 ; i<n ; i++){
        matriz_distancias[i].resize(n);
    }
  
   
    while(!ifile.eof()){
        
        ifile >> fila;
        ifile >> columna;
        ifile >> distancia;
        
        matriz_distancias[fila][columna] = distancia;
        matriz_distancias[columna][fila] = distancia;
        
    }

  /*  cout<< "Matriz : " << matriz_distancias.size();
    for(int i=0; i<n; i++){
       for(int j =0; j<n ; j++){
           cout << matriz_distancias[i][j] << " ";
       }
       cout << endl;
    }*/
}

void BL::generar_seleccionados(){
    Set_random(10);
    int x;
    vector<int> v;
    while(v.size() < m){
        x = Randint(0, n);
        if(std::find(v.begin(), v.end(), x) == v.end()) {
            v.push_back(x);
        } 
    }

    puntos p;
    
    for(int i=0; i<v.size(); i++){
        p.indice = v[i];
        sel.push_back(p);
       //cout << sel[i].indice << "   " << sel[i].contribucion <<endl;;
    }

    for(int i=0; i<n; i++){
        if(std::find(v.begin(), v.end(), i) == v.end()) {
            p.indice = i;
            vecinos.push_back(p);
        } 
        
    }

    std::random_shuffle(vecinos.begin(), vecinos.end());
  /*  vector<puntos>::iterator it;
    it = vecinos.begin();
    cout << "VECINOS NO EXPLORADOS: " <<endl;
    int contador =0;
    for(it; it != vecinos.end(); it++){
        cout << contador << " " << (*it).indice << endl;
        contador++;
    }*/

}

double BL::funcion_objetivo_total(vector<puntos> v){
    double suma=0;
    //cout<< "FILAS: " << matriz_distancias.size() << endl;
    //cout<< "COLUMNAS: " << matriz_distancias[0].size() <<endl;

    vector<puntos>::iterator it1, it2;
    it1 = it2 = v.begin();
    for(it1; it1 != v.end()-1; it1++){
        //cout << "entra   "<< m << endl;
        it2 = it1+1;
        for( it2; it2 != v.end(); it2++){
            suma += matriz_distancias[(*it1).indice][(*it2).indice];
        }
    }

    //cout << "Funcion objetivo: " << suma << endl;

    return suma;

}

vector<puntos> BL::ordenar_por_contribucion(vector<puntos> v){
    //Ordenamos de menor a mayor contibución

    puntos aux;
    for(int i=0; i<v.size(); i++){
        for(int j=0; j<v.size()-1; j++){
            if(v[j].contribucion > v[j+1].contribucion){
                aux = v[j];
                v[j] = v[j+1];
                v[j+1] = aux;
            }
        }
    }
    vector<puntos>::iterator it;
  /*  it = v.begin();
    for(it; it!=v.end(); it++){
       cout << (*it).indice << "   " << (*it).contribucion <<endl;;
    }*/

    return v;

}

void BL::generar_contribuciones(){
    for(int i=0; i< sel.size(); i++){
        sel[i].contribucion = funcion_objetivo(sel[i].indice, sel);
    }

}

void BL::factorizar(int anterior, int nuevo){
    vector<puntos>::iterator it;
    it = sel.begin();
    for(it; it!=sel.end(); it++){
        (*it).contribucion = (*it).contribucion - matriz_distancias[anterior][(*it).indice] + matriz_distancias[nuevo][(*it).indice] ;
    }

}
double BL::funcion_objetivo(int e, vector<puntos> v){
    double contribucion = 0;
    vector<puntos>::iterator it;
    it = v.begin();
    for(it; it != v.end(); it++){
        if((*it).indice != e){
            contribucion += matriz_distancias[(*it).indice][e];
        }
    }

    return contribucion;
}

puntos BL::genera_vecino(vector<puntos> v1, vector<puntos> v2){
    //Generamos un vecino que no se encuentre ni en v1 ni en v2
    //Lo usaremos para generar vecinos que no estén ni en seleccionados ni en explorados
    puntos vecino;
    vecino.indice =-1;
    puntos x;
    bool esta1 = false, esta2=false;
    vector<puntos>::iterator it1 , it2, it3;
    it1 = v1.begin();
    it2 = v2.begin();
    it3 = vecinos.begin();
    while(vecino.indice == -1 && it3 != vecinos.end()){
        x = (*it3);
        for(it1; it1!=v1.end(); it1++){
            if((*it1).indice == x.indice){
                esta1 = true;
            }
        }
        for(it2; it2!=v2.end(); it2++){
            if((*it2).indice == x.indice){
                esta2 = true;
            }
        }
        //cout << "Vecino " << vecino.indice << endl;
        if(!esta1 && !esta2){
            vecino = x;
            vecinos.erase(it3);
        }
        it3++;
      
    }

    

    return vecino;

}

void BL::setSeleccionados(vector<int> v){
    sel.clear();
    puntos p;
    vector<int>::iterator it;
    it = v.begin();
    for(it; it!=v.end(); it++){
        p.indice = (*it);
        sel.push_back(p);
    }


   /* vector<puntos>::iterator it2;
    it2 = sel.begin();
    cout << "SEL " <<endl;
    for(it2; it2!=sel.end(); it2++){
        cout << (*it2).indice << "  " ;
    }
    cout << "END" <<endl;*/

    for(int i=0; i<v.size(); i++){
        p.indice = v[i];
        sel.push_back(p);
    }

    vecinos.clear();

    for(int i=0; i<n; i++){
        if(std::find(v.begin(), v.end(), i) == v.end()) {
            p.indice = i;
            vecinos.push_back(p);
        } 
        
    }


    generar_contribuciones();
    sel = ordenar_por_contribucion(sel);

}

vector<int> BL::getSeleccionados_binario(){

    vector<int> indices_binario;

    for(int i=0; i<n; i++){
        indices_binario.push_back(0);
    }

    /*vector<puntos>::iterator it2;
    it2 = sel.begin();
    cout << "SEL " <<endl;
    for(it2; it2!=sel.end(); it2++){
        cout << (*it2).indice << "  " ;
    }
    cout << "END" <<endl;*/

    vector<puntos>::iterator it;
    it = sel.begin();
    int indice;
    for(it; it!=sel.end(); it++){
        indice = (*it).indice;
        indices_binario.at(indice) = 1;
    }

   /* vector<int>::iterator it3;
    it3 = indices_binario.begin();
    cout << "SEL BINARIO " <<endl;
    for(it3; it3!=indices_binario.end(); it3++){
        cout << (*it3) << "  " ;
    }
    cout << "END" <<endl;*/

    return indices_binario;
}