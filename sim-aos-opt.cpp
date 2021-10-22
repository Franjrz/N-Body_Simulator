#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <random>
#include <iomanip>

using namespace std;

struct objeto
{
  double p[3];
  double v[3];
  double m;
} ;

int fuerzaGravitatoria(double x1, double y1, double z1, double m1, double x2, double y2, double z2, double m2, double *resultado)
{
    //Si ambos objetos son los mismos entonces devolver 0
    if((x1==x2)&&(y1==y2)&&(z1==z2)&&(m1==m2))
        return -1;
    //Constante de gravitacion universal
    double G = 6.674e-11;
    //Calculo de la parte escalar (G*m1*m2)/||pj-pi||^3
    double pEscalar = (G*m1*m2)/pow(pow(x2-x1, 2.0) + pow(y2-y1, 2.0) + pow(z2-z1, 2.0), 3.0/2.0);
    //Calculo del vector resultado
    resultado[0] = (x2-x1)*pEscalar;
    resultado[1] = (y2-y1)*pEscalar;
    resultado[2] = (z2-z1)*pEscalar;
    return 0;
}

void vectorAceleracion(double x1, double y1, double z1, double m1, struct objeto *O, int sizeObjetos, double *resultado)
{
    //bucle for homologo al sumatorio de la formula
    for(int i = 0; i < sizeObjetos; i++){
        //Se calcula la fuerza que ejerce uno sobre el otro en forma vectorial
        double fuerza[3] = {0.0, 0.0, 0.0};
        int r = fuerzaGravitatoria(x1, y1, z1, m1, O[i].p[0], O[i].p[1], O[i].p[2], O[i].m, fuerza);
        //Se añade al resultado final, como indica el sumatorio
        if(r == 0){
            for(int j = 0; j < 3; j++){
                resultado[j] += fuerza[j];
            }
        }
    }
    //Se reescala por 1/m1
    double k = 1.0/m1;
    for(int j = 0; j < 3; j++){
        resultado[j] *= k;
    }
}

void colisiones(struct objeto *universo, int sizeObjetos)
{
    for (int j = 0; j < sizeObjetos; j++){
        if(j < sizeObjetos-1){
            for (int k = j+1; k < sizeObjetos; k++){
                if((universo[j].m > 0)&&(j != k)&&((pow(universo[k].p[0]-universo[j].p[0], 2.0) + pow(universo[k].p[1]-universo[j].p[1], 2.0) + pow(universo[k].p[2]-universo[j].p[2], 2.0)) < 1)){//Si 2 objetos son distintos y colisionan
                    universo[j].m += universo[k].m; //Se suman sus masas
                    universo[k].m = 0; //Se anula la masa del ultimo objeto
                    //Se suman las velocidades
                    for(int l = 0; l < 3; l++){
                        //Suma de velocidades
                        universo[j].v[l] += universo[k].v[l];
                    }
                }
            }
        }
    }
}

void guardarDatos(struct objeto *datos, string archivo, double size_enclosure, double time_step, int num_objects)
{
    ofstream Resultados(archivo);

    //Definicion de num decimales a 3
    Resultados << fixed << setprecision(3);

    //Cabecera
    Resultados << size_enclosure << " " << time_step << " " << num_objects << endl;

    //Cuerpo
    for (int i = 0; i < num_objects; i++){
        if(datos[i].m > 0){
            //Posiciones
            for(int j = 0; j < 3; j++){
                Resultados << datos[i].p[j] << " ";
            }
            //Velocidades
            for(int j = 0; j < 3; j++){
                Resultados << datos[i].v[j] << " ";
            }
            //Masa
            Resultados << datos[i].m << endl;
        }
    }

    //se cierra el archivo
    Resultados.close();
}

void simulation(const int num_objects, int num_iterations, int random_seed, double size_enclosure, double time_step)
{
    //Se inicializan los objetos de la simulation
    //random_device rd; //Se usa para obtener una semilla aleatoria con la que inicializar gen64
    mt19937_64 gen64(random_seed); //Generador de numeros pseudoaleatorios
    uniform_real_distribution<> disPos(0.0, size_enclosure); //Distribucion para generar las posiciones
    normal_distribution<> disMas{10e21, 10e15}; //Distribucion para generar las masas

    //elementos
    struct objeto universo[num_objects];

    //se generan de forma aleatoria los elementos
    for(int i = 0; i < num_objects; i++){
        objeto Obj;
        //posiciones y velocidades
        for(int j = 0; j < 3; j++){
            Obj.p[j] = disPos(gen64);//Se inicializan las posiciones de forma aleatoria
            Obj.v[j] = 0.0;//Se inicializan a 0 las velocidades
        }
        Obj.m = disMas(gen64)/10.0;//Se inicializan las masas de forma aleatoria
        universo[i] = Obj;
    }

    //Almacenar configuracion inicial
    guardarDatos(universo, "init_config.txt", size_enclosure, time_step, num_objects);

    //bucle de la simulation
    for(int i = 0; i < num_iterations; i++){
        cout << "Iteracion: " << i << endl;
        //comprobar si hay colisiones
        colisiones(universo, num_objects);
        //se crea la siguiente iteracion
        struct objeto newIt[num_objects];

        for (int j = 0; j < num_objects; j++){
            if(universo[j].m > 0){
                //Calculo del vector aceleracion
                double aceleracion[3] = {0.0, 0.0, 0.0};
                vectorAceleracion(universo[j].p[0], universo[j].p[1], universo[j].p[2], universo[j].m, universo, num_objects, aceleracion);
                //Calculo de las nuevas velocidades y posiciones en los 3 ejes
                objeto Obj;
                for(int k = 0; k < 3; k++){
                    Obj.v[k] = universo[j].v[k] + aceleracion[k]*time_step;
                    Obj.p[k] = universo[j].p[k] + Obj.v[k]*time_step;
                    //Calculo de rebotes
                    if(Obj.p[k] < 0.0){
                        Obj.p[k] = 0.0;
                        Obj.v[k] *= -1.0;
                    }
                    if(Obj.p[k] > size_enclosure){
                        Obj.p[k] = size_enclosure;
                        Obj.v[k] *= -1.0;
                    }
                    //cout << "rebotes calculados" << endl;
                }
                //cout << "Nueva masa de asteroide " << j << " " << universo.m[j] << endl;;
                //Cargar las masas en la nueva iteracion
                Obj.m = universo[j].m;
                newIt[j] = Obj;
            }
        }
        //cout << "Actualizando universo" << endl;
        //Actualizacion del universo
        universo = newIt;
        //cout << "Universo actualizado" << endl;
    }
    //Almacenar configuracion final
    guardarDatos(universo, "final_config.txt", size_enclosure, time_step, num_objects);
}

int main(int argc, char *argv[])
{
    cout << "sim-soa invoked with " << argc-1 << " parameters" << endl;

    string parametros[5] = {"num_objects", "num_iterations", "random_seed", "size_enclosure", "time_step"};

    for(int i = 0; i < 5; i++){
        cout << parametros[i] << ": ";
        if(i < argc-1)
            cout << argv[i+1] << endl;
        else
            cout << "?" << endl;
    }

    if(argc-1 != 5){
        cout << "Error: Numero de comandos incorrecto" << endl;
        exit(-1);
    }

    if(atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atof(argv[4]) <= 0.0 || atof(argv[5]) <= 0.0){
        cout << "Error: Algun argumento invalido" << endl;
        exit(-2);
    }
    const int num_objects = atoi(argv[1]);

    simulation(num_objects, atoi(argv[2]), atoi(argv[3]), atof(argv[4]), atof(argv[5]));

    return 0;
}
