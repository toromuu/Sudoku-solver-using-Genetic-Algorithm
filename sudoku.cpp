/* ------------------------- PROBLEMA DEL SUDOKU ----------------------- */
#include <ga/GASimpleGA.h>
#include <ga/GA1DArrayGenome.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

///PlantillaSudoku.txt
struct plantilla{
       int tam;
       int *fijo;
};

void leerSudoku(struct plantilla *S,char *nombreF){
   ifstream f(nombreF);

   f>>S->tam;

   S->fijo = new int[S->tam*S->tam];

   for(int i=0;i<S->tam*S->tam;i++)
           f>>S->fijo[i];

   f.close();
}

///MutacionSudoku.txt
bool checkColumna(int col[], int * check, int tam){
     bool repe=false;

     for(int i=0;i<tam;i++) check[i]=0;

     for(int i=0;i<tam;i++)
             check[col[i]-1]++;
     for(int i=0;i<tam;i++) if (check[i]>1) repe=true;

     return repe;
}

int MutacionSudoku(GAGenome& g,float pmut){

     // Hacer el casting correspondiente para obtener genome  //
     GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;


    struct plantilla * plantilla1;
    plantilla1 = (struct plantilla *) genome.userData();
    int nmut=0;
    int aux;
    //aux1;
    int fil;
    bool fila;

    int caux[plantilla1->tam];
    int *checkC=new int[plantilla1->tam];

    if (pmut<=0.0) return 0;

    for(int f=0; f<plantilla1->tam; f++)
       for(int c=0; c<plantilla1->tam; c++)
          if (plantilla1->fijo[(f*plantilla1->tam)+c]==0){
           if (GAFlipCoin(pmut) ){
                if (GAFlipCoin(0.5)) fila = true;
                else fila = false;

                if (!fila){

                      for(int j=0;j<plantilla1->tam;j++) caux[j]=genome.gene((j*plantilla1->tam)+c);
                      if (checkColumna(caux,checkC,plantilla1->tam)){
                         int v1 = GARandomInt(0,plantilla1->tam-1);
                         while (checkC[v1]<=1) v1=(v1+1)%plantilla1->tam;
                         v1++;
                         int v2 = GARandomInt(0,plantilla1->tam-1);
                         while (checkC[v2]!=0) v2=(v2+1)%plantilla1->tam;
                         v2++;

                         bool encontrado = false;
                         for(int j=0;j<plantilla1->tam && !encontrado;j++)
                                 if ((plantilla1->fijo[j*(plantilla1->tam)+c]==0)&&(genome.gene(j*(plantilla1->tam)+c)==v1)){
                                    encontrado = true;
                                    genome.gene((j*plantilla1->tam)+c,v2);
                                    fil = j;
                                 }

                         int col=(c+1)%plantilla1->tam;
                         while(genome.gene((fil*plantilla1->tam)+col)!=v2) col=(col+1)%plantilla1->tam;
                         if (plantilla1->fijo[(fil*plantilla1->tam)+col]==0) {
                                nmut++;
                                genome.gene((fil*plantilla1->tam)+col,v1);
                         }
                         else {
                              genome.gene((fil*plantilla1->tam)+c,v1);
                         }

                      }

                }
                else{
                   int v1 = (c + 1) %plantilla1->tam;
                   while ((plantilla1->fijo[(f*plantilla1->tam)+v1]!=0)) v1=(v1+1)%plantilla1->tam;
                   aux = genome.gene((f*plantilla1->tam)+c);
                   genome.gene((f*plantilla1->tam)+c,genome.gene((f*plantilla1->tam)+v1));
                   genome.gene((f*plantilla1->tam)+v1,aux);
                   nmut++;
                }
           }
          }

    return nmut;
}

///CruceSudoku.txt
int CruceSudoku(const GAGenome& p1,const GAGenome & p2,GAGenome* c1,GAGenome* c2){

     // Hacer el casting correspondiente para obtener m y p  //
     GA1DArrayAlleleGenome<int> & p = (GA1DArrayAlleleGenome<int> &)p1;
     GA1DArrayAlleleGenome<int> & m = (GA1DArrayAlleleGenome<int> &)p2;


    struct plantilla * plantilla1 = (struct plantilla *) m.userData();
    int n=0;

    int punto1=GARandomInt(0,m.length());
    while ((punto1%plantilla1->tam)!=0) punto1++;
    int punto2=m.length()-punto1;

    if (c1){
              // Hacer el casting correspondiente para obtener h1  //
              GA1DArrayGenome<int> & h1 = (GA1DArrayGenome<int> &)*c1;

             h1.copy(m,0,0,punto1); // el metodo copy esta definido en la clase GA1DArrayGenome
             h1.copy(p,punto1,punto1,punto2);
             n++;
    }

    if (c2){
              // Hacer el casting correspondiente para obtener h2  //
             GA1DArrayGenome<int> & h2 = (GA1DArrayGenome<int> &)*c2;

             h2.copy(p,0,0,punto1);
             h2.copy(m,punto1,punto1,punto2);
             n++;
    }

    return n;

}

///InicioSudoku.txt
void InicioSudoku(GAGenome& g){

     // Hacer el casting correspondiente para obtener genome  //
     GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;


     struct plantilla * plantilla1;
     plantilla1 = (struct plantilla *) genome.userData();

     int aux[plantilla1->tam];

     for(int f=0;f<plantilla1->tam;f++){

          for(int j=0;j<plantilla1->tam;j++) aux[j]=0;

          for(int j=1;j<=plantilla1->tam;j++){
            int v=GARandomInt(0,plantilla1->tam-1);
            while (aux[v]!=0) v=(v+1)%plantilla1->tam;
            aux[v]=j;
          }

          int i=0;

          while(i<plantilla1->tam){

              while((plantilla1->fijo[(f*plantilla1->tam)+i]==0) && (i<plantilla1->tam)) i++;

              if (i<plantilla1->tam){

                     bool encontrado=false;
                     for(int j=0;(j<plantilla1->tam) && (!encontrado);j++)
                             if (aux[j]==plantilla1->fijo[(f*plantilla1->tam)+i]) {
                                encontrado=true;
                                aux[j]=aux[i];
                             }

                     aux[i]=plantilla1->fijo[(f*plantilla1->tam)+i];
              }
              i++;

          }

          for(int c=0;c<plantilla1->tam;c++)
                  genome.gene((f*plantilla1->tam)+c,aux[c]);
     }
}

///Funcion de terminacion
GABoolean Termina(GAGeneticAlgorithm & ga){
    if ( (ga.statistics().minEver()==0) ||
        (ga.statistics().generation()==ga.nGenerations()) ) return gaTrue;
    else return gaFalse;
}


/* ------------------------- IMPLEMENTACION ----------------------- */

/* ---- Funcion Objetivo ----- */

/*

//// RECORRIDO POR EL GENOMA ////

	Cada individuo del sudoku se almacena como un array unidimensional, el genoma.

	La comprobación de las reglas del sudoku exigen la comprobación por filas,columnas y cajas.

	Para poder realizar estas comprobaciones hay que "cortar" el genoma segun el caso.

	El sudoku siempre será de un tamaño nxn, jugando con esto se pueden saber los limites inferiores
	y superiores de cada fila o columna.


	+ Para la división por filas.
		Con n*fila se puede determinar el limite superior de la fila actual y el inferior de la siguiente

		Ejemplo:  sudoku 3x3 cuyo genoma 123456789
		//se empieza siempre con fila=1

			1º fila	->  1 2 3  // 3*1 = 3
			2º fila ->	4 5 6  // 3*2 = 6
			3º fila ->	7 8 9  // 3*3 = 9


	+ Para la división por columnas
		Se va a aprovechar el hecho de que cada elemento siguiente del elemento que se esta tratando
		en la columna actual se encuentra a n posiciones de distancia en el array del genoma.
		Es decir para el calculo del valor de la posición i de la matriz seria
		valor_iesimo_columna j= genoma(i + n)


		Ejemplo:  sudoku 3x3 cuyo genoma 123456789
				  con disposición matricial seria:

						1 2 3
						4 5 6
						7 8 9
			Para obtener el valor "4", seria
			i = 0 + n --> genoma [3] = 4

			1º columna	->  1 4 7
			2º columna	->	2 5 8
			3º columna	->	3 6 9


	+ Para la división por cajas.
		El tablero de sudoku se puede dividir en cajas de tamaño igual  a raiz cuadrada de n
		Asi que si tuvieramos un tablero 4*4 , tendriamos 4 cajas cada una con 4 casillas
		En cada fila de cajas solo hay raiz cuadrada de n cajas, asi que cada vez que se llegue
		al final de la fila hay que bajar a la fila de cajas siguientes.


		1 2 | 3 4
		4 3 | 2 1   --> 1º fila de cajas
        ----|----
		3 1 | 4 2   --> 2º fila de cajas
		2 4 | 1 3

	    En cada caja hay que comprobar que no se repite ningun elemento.
		El recorrido dentro de la caja se hara de izquierda a derecha y de arriba a abajo.
		Para ello se aprovechara el tamaño de la caja, el cual servira de indicador para una
		vez recorrido los elementos de una fila de la caja, recorrer la fila de abajo.


		El calculo del i-elemento dentro del genoma se aprovecha del hecho de que los elementos
		de la siguiente fila dentro de una caja, se corresponden con los elementos en el genoma
		que se encuentran a un salto igual a  "tamaño de caja" de la posicion del ultimo elemento de la fila actual.

		Por ejemplo dado el siguiente genoma

		posicion 						 0 1 2 3 4 5 6 7
		genoma = 						 1 2 3 4 4 3 2 1
		tratamiento primera caja		 - -     - -

		Los elementos en la posicion i=0, i=1 se corresponden con los elementos de la primera fila
		de la primera caja, el primer elemento se encuentra a una distancia de +2 de la posicion del
		ultimo elemento de la actual fila, es decir
		posicion_primer_elemento_segunda_fila = 1 (de posicion actual) + 2 (de tamaño de caja) --> 3

		Por lo que el elemento que esta en la posicion 3 del genoma es el 4.




////CONTEO DE REPETCIONES////

	Para poder contar el número de repeticiones usaremos un array de n tamaño
	donde almacenaremos el número de apariciones de cada valor en su correspondiente
	i-elemento-1.

	Una vez obtenidas el numero de apariciones de la fila, columna, caja
	Se recorre el array de apariciones para buscar aquellos elementos que
	han aparecido mas de una vez, pues seran repeticiones.

	Se considera repeticion si el valor de la aparicion es por tanto > 1

*/

float Objective(GAGenome& g){



    //Definimos el genoma
    GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &)g;

	//Obtenemos n de la plantilla sudoku
    struct plantilla * plantillasudoku;
    plantillasudoku = (struct plantilla *) genome.userData();
    int n = plantillasudoku->tam;

    int apariciones[n] = { 0 }; // Para contar el número de repeticiones de cada valor
    int nrepeticiones = 0; // valor fitness a devolver en funcion del numero de elementos en apariciones con un valor >1

    int fila = 1;

	//genome.length() = n*n
	// El tratamiento por filas solo se tiene que hacer n veces
    for ( int i=0; i < genome.length(); i+=n ){

		// Se recorre secuencialmente la fila para contar
		// el número de veces que aparece un valor
		// este se almacena en el array de apariciones

        for ( int j = i; j < (n*fila); j++ ) {
			apariciones[genome.gene(j)-1]++;
		}

		// Se recorre el array de apariciones con los
		// valores actuales de la fila tratada
		// en busca de los elementos que tengan más de 1
		// aparicion, obteniendo el numero de repeticiones

        for (int z=0; z < n; z++){

            if ( (apariciones[z] - 1)  > 0 ){
				nrepeticiones += apariciones[z] - 1;
			}
            // Hay que limpiar el array de apariciones para la siguiente iteración
            apariciones[z] = 0;
        }

        fila++;
    }


    // El tratamiento por columnas es similar, solo que en este caso
	// Se va a aprovechar el hecho de que cada elemento siguiente del elemento que se esta tratando
	// en la columna actual se encuentra a n posiciones de distancia en el array del genoma.
    for( int i=0; i<n; i++ ){


        for( int j=i; j<genome.length() ; j+=n ){
            apariciones[genome.gene(j)-1]++;
		}

        for (int z=0; z < n; z++){

            if ( (apariciones[z] - 1)  > 0 ){
				nrepeticiones += apariciones[z] - 1;
			}
            // Hay que limpiar el array de apariciones para la siguiente iteración
            apariciones[z] = 0;
        }


    }


    //Para el tratamiento de cada caja
	// Obtenemos el tamano de cada caja con la raiz cuadrada de n
    int tamanoCaja = sqrt(n);

    int iCajaEnFilaActual = 0;

    for ( int posicionGenoma = 0; posicionGenoma < genome.length(); posicionGenoma += tamanoCaja ){

		// Cada fila
        if (iCajaEnFilaActual == tamanoCaja){
			//Hay que calcular la posicion del primer elemento de la primer fila
			// de la primera de caja de la fila de abajo
            posicionGenoma+=((tamanoCaja - 1)*n);
            iCajaEnFilaActual = 0;

			// Para evitar errores si al recalcular esta posicion se sobrepasa la longuitud del genoma
			// se devuelve el valor fitness
            if (posicionGenoma >= genome.length())
                return nrepeticiones;
        }

		// Calculamos el limite de posicion de la caja actual en el genoma
        int limitePosicionGenomaposicionGenoma = posicionGenoma + (n * (tamanoCaja - 1)) + tamanoCaja;


		// Nos sirve para determinar si se llego al ultimo elemento de la fila actual
        int limiteEFila = 0;

        for (int i = posicionGenoma; i < limitePosicionGenomaposicionGenoma; i++){

			//Si se llego al ultimo elemento de la fila actual en la iteracion pasada
            if ( limiteEFila == tamanoCaja){
				//Hay que calcular la posicion del primer elemento de la siguiente fila
                i += tamanoCaja * (tamanoCaja - 1);
				// Se vuelve a la "primera" columna
				limiteEFila = 0;
            }
			// Se accede al elemento i-esimo del genoma para recuperar el valor
            apariciones[genome.gene(i)-1]++;
            limiteEFila++;
        }

		// Contamos las repeticiones para esta caja
        for (int z=0; z < n; z++){

            if ( (apariciones[z] - 1)  > 0 ){
				nrepeticiones += apariciones[z] - 1;
			}
            // Hay que limpiar el array de apariciones para la siguiente iteración
            apariciones[z] = 0;
        }

        iCajaEnFilaActual++;
    }

    return nrepeticiones;
}

/* ---- Funcion Main ----- */
/*

*/
int main(int argc, char **argv){


// Creamos la plantila a partir del fichero de prueba dado
    char *ficheroPrueba = argv[1];
    struct plantilla *P = new struct plantilla;
    leerSudoku(P, ficheroPrueba);
// Obtenemos el tamano del caso
    int n = P->tam;

// Declaramos variables para los parametros del GA y las inicializamos
// a partir de los parametros de entrada

	int popsize = atoi(argv[2]);  // tamano de la poblacion {100,150}
    //int ngen = atoi(argv[3]); Para esta practica sera siempre 12000
	int ngen = 12000;             //numero generaciones
    float pcross = atof(argv[3]); // probabilidad de cruce {0.8, 0.85, 0.9, 0.95}
    float pmut = atof(argv[4]);   // probabilidad de mutacion {0.05, 0.075,0.1, 0.125}

	/*
    cout << "Parametros:    - Tamano poblacion: " << popsize << endl;
    cout << "               - Numero de generaciones: " << ngen << endl;
    cout << "               - Probabilidad cruce: " << pcross << endl;
    cout << "               - Probabilidad mutacion: " << pmut << endl << endl;
	*/

// Conjunto enumerado de alelos --> valores posibles de cada gen del genoma
// En el caso del sudoku desde 1-N, no se considera el 0
    GAAlleleSet<int> alelos;
    for(int i=1;i<n+1;i++) alelos.add(i);

// Creamos el genoma y definimos operadores de inicio, cruce y mutación
// (Tener en cuenta las estructuras y métodos dados)

    GA1DArrayAlleleGenome<int> genome(n*n,alelos,Objective,P); //indicar la plantilla P

	genome.crossover(CruceSudoku);
    genome.mutator(MutacionSudoku);
	genome.initializer(InicioSudoku);

// Creamos el algoritmo genetico
    GASimpleGA ga(genome);

// Inicializamos - minimizar funcion objetivo, tamaño poblacion, nº generaciones,
// pr. cruce y pr. mutacion, selección y le indicamos que evolucione.

    ga.minimaxi(-1);
    ga.populationSize(popsize);
    ga.nGenerations(ngen);
    ga.pCrossover(pcross);
    ga.pMutation(pmut);

	// Declaramos los selectores posibles en función del parametro indicado
	// 1 ruleta , 2 torneo
	string selector;

	if (atoi(argv[5])== 1){
		GARouletteWheelSelector selectorRuleta;
		ga.selector(selectorRuleta);
		selector = "Ruleta";// Para la impresion del operador utilizado en el fichero
	}
	else if (atoi(argv[5])== 2) {
			GATournamentSelector selectorTorneo;
			ga.selector(selectorTorneo);
			selector = "Torneo";
		}

    ga.terminator(Termina);

	//obligatoriamente utilizaremos en valor 1 en el método (evolve(1)).
    ga.evolve(1);

/* ---- Exportacion de los resultados ----- */

//  Por un lado crearemos un fichero salida  por cada fichero-prueba distinto
//  Cada combinacion de parametros se anadira al final del fichero, por lo que para evitar datos duplicados para el analisis
//  es importante borrar el fichero antes de ejecutar el script.

    std::ofstream ficheroSalida;
    char ficheroSalida1[50] = "ficheroSalida-";
    strcat(ficheroSalida1, ficheroPrueba); //Concatenacion para id unico de fichero

    ficheroSalida.open(ficheroSalida1,std::ofstream::in | std::ofstream::out | std::ofstream::app);

 // Volcamos en el fichero de salida los parametros utilizados en la ejecucion
    ficheroSalida << "Fichero: " << ficheroPrueba << endl;
	ficheroSalida << "Parametros: "<< endl;
    ficheroSalida << "+ Tamano poblacion: " << popsize << endl;
    ficheroSalida << "+ Numero de generaciones: " << ngen << endl;
    ficheroSalida << "+ Probabilidad cruce: " << pcross << endl;
    ficheroSalida << "+ Probabilidad mutacion: " << pmut << endl;
    ficheroSalida << "+ Selector: " << selector << endl << endl;


// Escribimos en el fichero salida el resultado obtenido
    ficheroSalida << "El AG ha encontrado la siguiente solucion " << endl;

    GA1DArrayGenome<int> & bestGenome = (GA1DArrayGenome<int> &) ga.statistics().bestIndividual();

	for (int i = 0; i < bestGenome.length(); i+=n){
        for (int j = i; j < i+n; j++){
            if (j == i+n-1)
                ficheroSalida << bestGenome.gene(j);
            else
                ficheroSalida << bestGenome.gene(j) << " ";
        }
        ficheroSalida << endl;
    }

    ficheroSalida << endl;
    ficheroSalida << "valor fitness " << ga.statistics().minEver() << endl << endl << endl ;

 // Por otro lado para el analisis de los datos y la generacion de la tabla resultado con todos los casos posibles
 // Crearemos un documento excel
     std::ofstream ficheroCsv;
    char excel[50] = "tablaResultado-";
    strcat(excel, ficheroPrueba);
    strcat(excel, ".csv");
    ficheroCsv.open(excel,std::ofstream::in | std::ofstream::out | std::ofstream::app);

 // Volcamos en el fichero de salida los parametros utilizados en la ejecucion
    //ficheroCsv << "Selector"<<","<<"TPoblacion"<<","<<"PCruce"<<","<<"PMutacion"<<","<<"Fitness"<< endl;
    ficheroCsv << selector<<","<<popsize<<","<<pcross<<","<<pmut<<","<<ga.statistics().minEver()<< endl;


}


