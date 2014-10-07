#include <stdio.h>
#include <stdlib.h>

int main(){

  FILE *nombres = fopen("nombres_datos.dat", "r");
  FILE *nombresmodificados = fopen("nombres_modificados.sh", "w");
  char titulo;
  int n = 0;
  
  while ((titulo = getc(nombres)) != EOF){
    ungetc(titulo, nombres);
    n = 0;
    char datos[100] = "";
    while ((titulo = getc(nombres)) != '\n' && titulo != EOF){
      datos[n] = titulo;
      n++;
    }
    fprintf(nombresmodificados, "%s%s \n", "./a.out ", datos);
  }
  fclose(nombres);
  fclose(nombresmodificados);
  return 0;
}
