# Calculo Numerico Projeto UFPE
Projeto utilizando os métodos numéricos: Bisseção, Newton-Raphson e Halley

## Exemplo de uso da biblioteca


```C++
     int main(){
        int metodo, Nmax, i;
        double e1, e2, c[6], a, b, media;
        
        scanf("%d", &metodo); // escaneando o método

        for(i=0; i<6; i++){  //constantes 
          scanf("%lf", &c[i]);
        }

        scanf("%d %lf %lf %lf %lf", &Nmax, &e1, &e2, &a, &b); //resto das entradas

        media = (a+b)/2;

        if(metodo == 1){     //a escolha do método chama a respectiva função
          bissecao(a, b, Nmax, e1, e2, c);
        }
        else if(metodo == 2){
          newton(media, Nmax, e1, e2, c);
        }
        else if(metodo == 3){
          halley(media, Nmax, e1, e2, c);
        }

        return 0;
      }
```
