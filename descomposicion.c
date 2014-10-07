#include<stdio.h>
#include<stdlib.h>
#include<math.h>


float *llenar_matriz(char *archivo, int numero_datos);
void llenar_coeficientes(float *matriz_datos, float *coeficientes, float *vectorb, float *tiempo, int n_row);
void eliminar_ruido(float *b);
float *transpuesta(float *coeficientes, int n_row, int n_cols);
float *multiplicacion(float *matrizA, float *matrizB, int n_rowA, int n_colsA, int n_rowB, int n_colsB);
void descomposicion_lu(float *matriz, float *b, float *U, float *L, int n);
void solucion_u(float *matriz, float *vectorb, float *y_0, float *v_y, float *g_y);
void imprimir_resultados(char *archivo, float y_0, float v_y, float g_y);

int main(int argc, char **argv){
  int n_cols = 3;
  char *archivo = argv[1];
  int n_row = 38;
  float *matriz_datos;
  matriz_datos = llenar_matriz(archivo, n_row);
  
  float *coeficientes = malloc(n_row*n_cols*sizeof(float));
  float *vectorb = malloc(n_row*sizeof(float));
  float *tiempo = malloc(n_row*sizeof(float));
  llenar_coeficientes(matriz_datos, coeficientes, vectorb, tiempo, n_row);
  eliminar_ruido(vectorb); 
  
  float *transpuesta_coeficientes;
  transpuesta_coeficientes = transpuesta(coeficientes,n_row,n_cols);
  
  float *matriz;
  float *b_prima;
  float *U = malloc(n_cols*n_cols*sizeof(float));
  float *L = malloc(n_cols*n_cols*sizeof(float));
  float y_0,v_y,g_y;
  
  matriz = multiplicacion(transpuesta_coeficientes, coeficientes, n_cols, n_row, n_row, n_cols);
  b_prima = multiplicacion(transpuesta_coeficientes, vectorb, n_cols, n_row, n_row, 1);
  descomposicion_lu(matriz, b_prima, U, L, n_cols);
  solucion_u(U, b_prima, &y_0, &v_y, &g_y);
  imprimir_resultados(archivo, y_0, v_y, g_y);
  return 0;
}

// Cada archivo tiene 38 datos de tiempo y posicion vertical. En esta funcion se toma en cuenta que el numero de datos es el numero de columnas de la matriz
float *llenar_matriz(char *archivo, int numero_datos){
  
  float *matriz;
  FILE *arch;
  int n_row = numero_datos, n_cols = 2;
  int i, j;

  if(!(arch = fopen(archivo, "r"))){
    printf("Problema abriendo el archivo %s\n", archivo);
    exit(1);
  }

  matriz = malloc(n_row*n_cols*sizeof(float));

  for (i = 0; i < n_row; i++){
    for (j = 0; j < n_cols; j++){
      fscanf(arch, "%f", &matriz[i*n_cols + j]);
    }
  }
  fclose(arch);
  return matriz;
}
      
// Un modelo que puede explicar un sistema fisico, se representa por una matriz de coeficientes que multiplica los parametros del modelo. En este caso se usa la ecuacion de movimiento parabolico en y, para conocer los coeficientes que multiplican los parametros y_0, v_y y g_y, que son funcion del tiempo.

void llenar_coeficientes(float *matriz_datos, float *coeficientes, float *vectorb, float *tiempo, int n_row){

  int i, j;
  for (i = 0; i < n_row; i++){
    coeficientes[i * 3] = 1;
    tiempo[i] = matriz_datos[i * 2];
    coeficientes[i * 3 + 1] = tiempo[i];
    coeficientes[i * 3 + 2] = -pow(tiempo[i],2)/2.0;
    vectorb[i] = matriz_datos[i * 2 + 1];
  }
}

// El ruido se caracteriza por tener en algun momento una sucesion de minimos y maximos locales, lo que no ocurre en una curva parabolica. Se eliminan los datos que tienen mas de tres maximos seguidos.

void eliminar_ruido(float *b){

  int i, maximos=0;
  for(i = 1; i < 37; i++){
    if (b[i] < b[i+1] && b[i] < b[i-1]){
      maximos++;
    }
    if(maximos > 3){
      exit(1);
    }
  }
}

// Para encontrar la mejor estimacion de los parametros del modelo, se hace uso de la transpuesta de los coeficientes encontrados con la funcion previa.

float *transpuesta(float *coeficientes, int n_row, int n_cols){

  float *t = malloc(n_row * n_cols * sizeof(float));
  int i, j;

  for (i = 0; i < n_row; i++){
    for (j = 0; j < n_cols; j++){
      t[j * n_row + i] = coeficientes[i * n_cols + j];
    }
  }
  return t;
}

// Para el metodo de descomposicion se usa extensivamente la multiplicacion de matrices.

float *multiplicacion(float *matrizA, float *matrizB, int n_rowA, int n_colsA, int n_rowB, int n_colsB){

  if (n_colsA != n_rowB){
    printf("NO se pueden multiplicar estas matrices");
    exit(1);
  }

  float *matrizM;
  int i, j, k, filasM = n_rowA, columnasM = n_colsB;
  float multiplicacion;
  
  matrizM = malloc (filasM * columnasM * sizeof(float));
  for (i = 0; i < filasM; i++){
    for (j = 0; j < columnasM; j++){
      multiplicacion = 0.0;
      for (k = 0; k < n_colsA; k++){
	multiplicacion += matrizA[i * n_colsA + k] * matrizB[k * n_colsB + j];
      }
      matrizM[columnasM * i + j] = multiplicacion;
    }
  }
  return matrizM;
}

// A continuacion se muestra el codigo fuente de la descomposicion de Cholesky, teniendo en cuenta que lo primero que se hace es revisar si se requiera cambiar las filas y si es el caso, las cambia.

void descomposicion_lu(float *matriz, float *b, float *U, float *L, int n){

  int i, j, k;
  int cont = 0;
  float fila = 0.0;

  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      U[i * n + j] = matriz[i * n + j];
    }
  }
  
  for (i = 0; i < n; i++){
    L[i * n + i] = 1;
    if (U[i * n + i] < 1.0E-10){
      for (j = i + 1; j < n && cont == 0; j++){
	if (U[j * n + i] != 0){
	  U[i * n + i] = U[j * n + i];
	  fila = b[i];
	  b[i] = b[j];
	  b[j] = fila;
	  cont = 1;
	  for (k = i + 1; k < n; k++){
	    fila = U[i * n + k];
	    U[i * n + k] = U[j * n + k];
	    U[j * n + k] = fila;
	  }
	}
      }
    }

    for (j = i + 1; j < n; j++){
      L[j * n + i] = U[j * n + i]/U[i * n + i];
      U[j * n + i] = 0;
      b[j] = b[j] - L[j * n + i] * b[i];
      for (k = i + 1; k < n; k++){
	U[j * n + k] = U[j * n + k] - L[j * n + i] * U[i * n + k];
      }
    }
  }
}


// Solucion del sistema de ecuaciones a partir de una matriz superior triangular resultado de la descomposicion de Cholesky.

void solucion_u(float *matriz, float *vectorb, float *y_0, float *v_y, float *g_y){

  int n = 3;
  int i,j;
  float producto;
  float *s = malloc(n * sizeof(float));
  s[n - 1] = vectorb[n - 1]/matriz[n*n - 1];

  for (i = n - 2; i > -1; i--){
    producto = 0.0;
    for (j = n - 1; j > i; j--){
      producto += s[j] * matriz[n * i + j];
    }
    s[i] = (vectorb[i] - producto) / (matriz[n * i + i]);
  }
  *y_0 = s[0];
  *v_y = s[1];
  *g_y = s[2];
}

// Imprimir los resultados en el archivo de salida, teniendo en cuenta encontrar los angulos phi y theta del nombre del archivo de datos

void imprimir_resultados(char *archivo, float y_0, float v_y, float g_y){

  int i = 0, j, angulo_phi = 0, angulo_theta = 0, unidad;
  float phi = 0.0, theta = 0.0;
  char *archivo_final = "resultados.dat";
  FILE *resultados = fopen(archivo_final, "a");

  while (archivo[i] != '\0' && (angulo_phi == 0 || angulo_theta == 0)){
    i += 1;
    if (angulo_theta == 0 && archivo[i] == 'e' && archivo[i + 1] == 't' && archivo[i + 2] == 'a'){
      angulo_theta = 1;
      j = i + 4;
      while (archivo[j] != '.'){
	unidad = archivo[j] - '0';
	theta = theta * 10.0 + unidad;
	j += 1;
      }
      j += 1;
      while (archivo[j] != '_'){
	unidad = archivo[j] - '0';
	theta = theta + unidad / 10.0;
	j += 1;
      }
    }
    if (angulo_phi == 0 && archivo[i] == 'p' && archivo[i + 1] == 'h' && archivo[i + 2] == 'i'){
      angulo_phi = 1;
      j = i + 4;
      while (archivo[j] != '.'){
	unidad = archivo[j] - '0';
	phi = phi * 10.0 + unidad;
	j += 1;
      }
      j += 1;
      while (archivo[j] != '.'){
	unidad = archivo[j] - '0';
	phi = phi + unidad / 10.0;
	j += 1;
      }
    }
  }
  fprintf(resultados, "%f %f %f %f %f \n", theta, phi, y_0, v_y, g_y);
  fclose(resultados);
}
