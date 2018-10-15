#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double resFuncao(double x, double* c);
double resDerivada(double x, double* c);
double resDerivadaSegunda(double x, double* c);
void bissecao(double a, double b, int Nmax, double e1, double e2, double* c);
void newton(double x, int Nmax, double e1, double e2, double* c);
void halley(double x, int Nmax, double e1, double e2, double* c);

double resFuncao(double x, double* c){

  return c[0] * cos(c[1] * x) + c[2] * sin(c[3] * x) + exp(c[4] * x) + c[5];
}

double resDerivada(double x, double* c){

  return -c[0] * sin(c[1] * x)*c[1] + c[2] * cos(c[3] * x)*c[3] + exp(c[4] * x)*c[4];
}

double resDerivadaSegunda(double x, double* c){
  
  return -c[0] * cos(c[1]*x) * c[1] * c[1] - c[2] * sin(c[3] * x) * c[3] * c[3] + exp(c[4] * x) * c[4] * c[4];
}

void bissecao(double a, double b, int Nmax, double e1, double e2, double* c){
  int i;
  double res1, res2, res, x;
  res1 = resFuncao(a, c);
  res2 = resFuncao(b, c);
  
  if(res1 > 0 && res2 < 0){    // compara qual dos valores é positivos e negativos
    for(i = 1; i <= Nmax; i++){
      x = (a+b)/2;
      
      res = resFuncao(x, c);

      
      if(res > 0){
        a = x;
      }
      if(res < 0){
        b = x;
      }

      
      if(b-a < e1 && fabs(res) < e2){
        printf("Metodo: Bissecao\n");
        printf("Iteracao: %d\n", i);
        printf("CONVERGIU\n");
        printf("Raiz final:%0.10lf\n", x);
        printf("|a-b|:%0.10lf\n", b-a);
        printf("|f(x_i)|:%0.10lf", fabs(res));
      break;
      }
      printf("Metodo: Bissecao\n");
      printf("Iteracao: %d\n", i);
      printf("NAO CONVERGIU\n");
      printf("Raiz final:%0.10lf\n", x);
      printf("|a-b|:%0.10lf\n", b-a);
      printf("|f(x_i)|:%0.10lf\n\n", fabs(res));
    }  
  }
  else if(res1 < 0 && res2 > 0){
    for(i = 1; i <= Nmax; i++){
      x = (a+b)/2;
      res = resFuncao(x, c);
      
      if(res > 0){
        b = x;
      }
      if(res < 0){
        a = x;
      }
        if(b-a < e1 && fabs(res) < e2){
          printf("Metodo: Bissecao\n");
          printf("Iteracao: %d\n", i);
          printf("CONVERGIU\n");
          printf("Raiz final:%0.10lf\n", x);
          printf("|a-b|:%0.10lf\n", b-a);
          printf("|f(x_i)|:%0.10lf", fabs(res));
          break;
      }

      printf("Metodo: Bissecao\n");
      printf("Iteracao: %d\n", i);
      printf("NAO CONVERGIU\n");
      printf("Raiz final:%0.10lf\n", x);
      printf("|a-b|:%0.10lf\n", b-a);
      printf("|f(x_i)|:%0.10lf\n\n", fabs(res));
    }  
  }
  else{  //se os dois valores tiverem o mesmo sinal, não dá para usar bisseção
    printf("Função e intervalo inválidos para o método\n");
    exit(1);
  }
}

void newton(double x, int Nmax, double e1, double e2, double* c){
  int i;
  double ab, res;
  for(i = 1; i <= Nmax; i++){
    ab = fabs(resFuncao(x, c)/resDerivada(x, c));
    x = x - resFuncao(x, c)/resDerivada(x, c);
    res = resFuncao(x, c);
    if(fabs(res) < e2 && fabs(ab) < e1){            
      printf("Metodo: Newton-Raphson\n");
      printf("Iteracao: %d\n", i);
      printf("CONVERGIU\n");
      printf("Raiz final:%0.10lf\n", x);
      printf("|x_i-x_(i-1)|:%0.10lf\n", ab);
      printf("|f(x_i)|:%0.10lf", res);
      break;
    }
    printf("Metodo: Newton-Raphson\n");
    printf("Iteracao: %d\n", i);
    printf("NAO CONVERGIU\n");
    printf("Raiz final:%0.10lf\n", x);
    printf("|x_i-x_(i-1)|:%0.10lf\n", ab);
    printf("|f(x_i)|:%0.10lf\n\n", res);
  } 
}

void halley(double x, int Nmax, double e1, double e2, double* c){
  int i;
  double ab, res;
  for(i = 1; i <= Nmax; i++){
    ab = fabs((2*resFuncao(x, c)*resDerivada(x, c))/(2 * pow(resDerivada(x, c), 2) - resFuncao(x, c) * resDerivadaSegunda(x, c)));
    x = x - (2*resFuncao(x, c)*resDerivada(x, c))/(2 * pow(resDerivada(x, c), 2) - resFuncao(x, c) * resDerivadaSegunda(x, c));
    res = resFuncao(x, c);
    if(fabs(res) < e2 && fabs(ab) < e1){            
          printf("Metodo: Halley\n");
          printf("Iteracao: %d\n", i);
          printf("CONVERGIU\n");
          printf("Raiz final:%0.10lf\n", x);
          printf("|x_i-x_(i-1)|:%0.10lf\n", ab);
          printf("|f(x_i)|:%0.10lf", res);
      break;
    }
    printf("Metodo: Halley\n");
    printf("Iteracao: %d\n", i);
    printf("NAO CONVERGIU\n");
    printf("Raiz final:%0.10lf\n", x);
    printf("|x_i-x_(i-1)|:%0.10lf\n", ab);
    printf("|f(x_i)|:%0.10lf\n\n", res);
  }   
}