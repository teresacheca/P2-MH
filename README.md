# P2-MH

Práctica 2: Maximum diversity problem (MDP)

Para el desarrollo de la práctica se ha creado un fichero distinto para cada variante del algoritmo. Por lo que
tendríamos un fichero para AGE, otro para AGG, otro para AGM y otro para BL con sus respectivos
ficheros .h.

Por lo que encontraríamos:

  -Estacionario:

    -estacionario.h: fichero donde están declarados todas las funciones para el algoritmo AGE
    -estacionario.cpp: fichero donde están implementadas todas las funciones del algoritmo AGE

-Generacional:
    
    -generacional.h: fichero donde están declarados todas las funciones para el algoritmo AGG
    -generacional.cpp: fichero donde están implementadas todas las funciones del algoritmo AGG

-Memético:

    -memetico.h: fichero donde están declarados todas las funciones para el algoritmo AGM
    -memético.cpp: fichero donde están implementadas todas las funciones del algoritmo AGM

Dependiendo de los parámetros de entrada, los algoritmos genéricos utilizarán el operador de cruce
uniforme o el operador de cruce basado en posición.

De la misma forma, dependiendo de los parámetros de entrada, los algoritmos Meméticos pueden
ser:

    AM-(10,1.0): Cada 10 generaciones, aplicar la BL sobre todos los cromosomas de la población
    AM-(10,0.1): Cada 10 generaciones, aplicar la BL sobre un subconjunto de cromosomas de la población s
    eleccionado aleatoriamente con probabilidad pLS igual a 0.1 para cada cromosoma
    AM-(10,0.1mej): Cada 10 generaciones, aplicar la BL sobre los 0.1·N mejores cromosomas de la
    población actual (N es el tamaño de ésta)

-Búsqueda Local:

    -BL.h: fichero donde están declarados todas las funciones para el algoritmo BL
    -BL.cpp: fichero donde están implementadas todas las funciones del algoritmo BL

Luego, encontraremos también el fichero Makefile y un script llamado ejecutable.
Este script hace una llamada por cada fichero de datos que se nos han proporcionado.
De esta forma, no tendremos nada más que ejecutar el script y de forma automática se
hará una ejecución del greedy y del algoritmo de búsqueda local para cada fichero de
datos.
estos ficheros de datos, es necesario que estén fuera del directorio software y dentro
un directorio llamado ‘Instancias y Tablas MDP 2020-21’

Por otra parte, la ejecución del programa es muy sencilla, para ello sólo hay que seguir
los siguientes pasos:

    1. Abrir una terminal y situarnos dentro de la carpeta software
    2. Escribir el comando: make
    3. Ejecutar el script ejecutable con el comando: ./ejecutable

Una vez hecho esto, se nos mostrarán por pantalla los resultados de los algoritmos
BL, AGG-uniforme, AGG-posición, AGE-uniforme, AGE-posición, AM-(10,1.0), AM-(10,0.1), AM-
(10,0.1mej), así como el tiempo que tardan en ejecutarse y la desviación típica de cada uno.
