#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <random>
#include <iomanip>

using namespace std;

struct objetos
{
  vector<double> p[3];
  vector<double> v[3];
  vector<double> m;
} ;

int fuerzaGravitatoria(double x1, double y1, double z1, double m1, double x2, double y2, double z2, double m2, vector<double> &resultado)
{
    //Si ambos objetos son los mismos entonces devolver 0
    if((x1==x2)&&(y1==y2)&&(z1==z2)&&(m1==m2))
        return -1;
    //Constante de gravitacion universal
    double G = 6.674e-11;
    //Calculo de la parte escalar (G*m1*m2)/||pj-pi||^3
    double pEscalar = (G*m1*m2)/pow(pow(x2-x1, 2.0) + pow(y2-y1, 2.0) + pow(z2-z1, 2.0), 3.0/2.0);
    //Calculo del vector resultado
    //cout <<  endl;
    //cout << "(" << x1 <<  ", " << y1 <<  ", " << z1 <<  ", " << m1 <<  ") (" << x2 <<  ", " << y2 <<  ", " << z2 <<  ", " << m2 <<  ")" << endl;
    //cout << pow(pow(x2-x1, 2.0) + pow(y2-y1, 2.0) + pow(z2-z1, 2.0), 1.0/2.0) << endl;
    resultado.push_back((x2-x1)*pEscalar);
    //cout << "Fuerza x " << (x2-x1)*pEscalar << endl;
    resultado.push_back((y2-y1)*pEscalar);
    //cout << "Fuerza y " << (y2-y1)*pEscalar << endl;
    resultado.push_back((z2-z1)*pEscalar);
    //cout << "Fuerza z " << (z2-z1)*pEscalar << endl;
    //cout <<  endl;
    return 0;
}

void vectorAceleracion(double x1, double y1, double z1, double m1, struct objetos &O, vector<double> &resultado, long unsigned int pos, ofstream& Fuerzas, ofstream& Aceleraciones)
{
    //bucle for homologo al sumatorio de la formula
    for(long unsigned int i = 0; i < O.p[0].size(); i++){
        //Se calcula la fuerza que ejerce uno sobre el otro en forma vectorial
        vector<double> fuerza;
        int r = fuerzaGravitatoria(x1, y1, z1, m1, O.p[0][i], O.p[1][i], O.p[2][i], O.m[i], fuerza);
        //Se añade al resultado final, como indica el sumatorio
        if(r == 0){
            transform(resultado.begin(), resultado.end(), fuerza.begin(), resultado.begin(), plus<double>());
            //cout << "fuerza (" << fuerza[0] <<  ", " << fuerza[1] <<  ", " << fuerza[2] <<  "), resultado (" << resultado[0] <<  ", " << resultado[1] <<  ", " << resultado[2] <<  ")" << endl;
        }
    }
    Fuerzas << "Force(" << pos << ") = " << resultado[0] << " " << resultado[1] << " " << resultado[2] << endl;
    //Se reescala por 1/m1
    double k = 1.0/m1;
    transform(resultado.begin(), resultado.end(), resultado.begin(), [k](double &c){ return c*k; });
    Aceleraciones << "Acceleration(" << pos << ") = " << resultado[0] << " " << resultado[1] << " " << resultado[2] << endl;
}

void colisiones(struct objetos &universo, ofstream& Posiciones)
{
    Posiciones << "DENTRO DE COLISIONES " << endl;
    for (long unsigned int j = 0; j < universo.p[0].size(); j++){
        if(j < universo.p[0].size()-1){
            for (long unsigned int k = j+1; k < universo.p[0].size(); k++){
                if((j != k)&&((pow(universo.p[0][k]-universo.p[0][j], 2.0) + pow(universo.p[1][k]-universo.p[1][j], 2.0) + pow(universo.p[2][k]-universo.p[2][j], 2.0)) < 1)){//Si 2 objetos son distintos y colisionan
                    Posiciones << "Collapsing objects " << j << " and " << k << endl;
                    Posiciones << "Body " << j << ": " << universo.p[0][j] << " " << universo.p[1][j] << " " << universo.p[2][j] << " " << universo.v[0][j] << " " << universo.v[1][j] << " " << universo.v[2][j] << " " << universo.m[j] << endl;
                    Posiciones << "Body " << k << ": " << universo.p[0][k] << " " << universo.p[1][k] << " " << universo.p[2][k] << " " << universo.v[0][k] << " " << universo.v[1][k] << " " << universo.v[2][k] << " " << universo.m[k] << endl;
                    //cout << "Se va a eliminar asteroide " << k << endl;
                    universo.m[j] += universo.m[k]; //Se suman sus masas
                    //cout << "Masa nueva de asteroide " << j << " " << universo.m[j] << endl;
                    //Se suman las velocidades y se destruye el segundo objeto
                    for(int l = 0; l < 3; l++){
                        //Suma de velocidades
                        //cout << "Velocidad " << l << " de " << j << " " << universo.v[l][j] << endl;
                        //cout << "Velocidad " << l << " de " << k << " " << universo.v[l][k] << endl;
                        universo.v[l][j] += universo.v[l][k];
                        //cout << "Velocidad " << l << " nueva de asteroide " << j << " " << universo.v[l][j] << endl;
                        //destruccion del segundo objeto
                        universo.p[l].erase(universo.p[l].begin()+k);
                        universo.v[l].erase(universo.v[l].begin()+k);
                    }
                    universo.m.erase(universo.m.begin()+k);
                    Posiciones << "Body " << k << " removed" << endl;
                    Posiciones << "Body " << j << " after collapse" << endl;
                    Posiciones << "Body " << j << ": " << universo.p[0][j] << " " << universo.p[1][j] << " " << universo.p[2][j] << " " << universo.v[0][j] << " " << universo.v[1][j] << " " << universo.v[2][j] << " " << universo.m[j] << endl;
                    Posiciones << endl;
                    //Se echa hacia atrás el contador k para no saltarse un asteroide
                    k --;
                    //cout << "Numero total de asteroides " << universo.p[0].size() << endl;
                }
            }
        }
    }
}

void guardarDatos(struct objetos &datos, string archivo, double size_enclosure, double time_step, int num_objects)
{
    ofstream Resultados(archivo);

    //Definicion de num decimales a 3
    Resultados << fixed << setprecision(3);

    //Cabecera
    Resultados << size_enclosure << " " << time_step << " " << num_objects << endl;

    //Cuerpo
    for (long unsigned int i = 0; i < datos.p[0].size(); i++){
        //Posiciones
        for(int j = 0; j < 3; j++){
            Resultados << datos.p[j][i] << " ";
        }
        //Velocidades
        for(int j = 0; j < 3; j++){
            Resultados << datos.v[j][i] << " ";
        }
        //Masa
        Resultados << datos.m[i] << endl;
    }

    //se cierra el archivo
    Resultados.close();
}

void simulation(int num_objects, int num_iterations, int random_seed, double size_enclosure, double time_step)
{
    //Se inicializan los objetos de la simulation
    //random_device rd; //Se usa para obtener una semilla aleatoria con la que inicializar gen64
    mt19937_64 gen64(random_seed); //Generador de numeros pseudoaleatorios
    uniform_real_distribution<> disPos(0.0, size_enclosure); //Distribucion para generar las posiciones
    normal_distribution<> disMas{10e21, 10e15}; //Distribucion para generar las masas

    //elementos
    struct objetos universo;

    //se generan de forma aleatoria los elementos
    for(int i = 0; i < num_objects; i++){
        //posiciones y velocidades
        for(int j = 0; j < 3; j++){
            universo.p[j].push_back(disPos(gen64));//Se inicializan las posiciones de forma aleatoria
            universo.v[j].push_back(0.0);//Se inicializan a 0 las velocidades
        }
        universo.m.push_back(disMas(gen64)/10.0);//Se inicializan las masas de forma aleatoria
    }

    ofstream Fuerzas("forces.txt");
    ofstream Aceleraciones("accelerations.txt");
    ofstream Velocidades("velocities.txt");
    ofstream Posiciones("positions.txt");

    for (long unsigned int j = 0; j < universo.p[0].size(); j++){
        Posiciones << "Position(" << j << ") = " << universo.p[0][j] << " " << universo.p[1][j] << " " << universo.p[2][j] << endl;
    }


    //Almacenar configuracion inicial
    guardarDatos(universo, "init_config.txt", size_enclosure, time_step, num_objects);

    //bucle de la simulation
    for(int i = 0; i < num_iterations; i++){

        Fuerzas << "\nIteration: " << i << "\nNet Forces" << endl;
        Aceleraciones << "\nIteration: " << i << endl;
        Velocidades << "\nIteration: " << i << endl;
        Posiciones << "\nIteration: " << i << endl;

        //comprobar si hay colisiones
        colisiones(universo, Posiciones);

        //se crea la siguiente iteracion
        struct objetos newIt;

        //cout << "Creado nuevo universo" << endl;

        for (long unsigned int j = 0; j < universo.p[0].size(); j++){
            //cout << "Calculando vector de aceleracion de asteroide " << j << endl;
            //Calculo del vector aceleracion
            vector<double> aceleracion(3, 0.0);
            //cout <<  endl;
            vectorAceleracion(universo.p[0][j], universo.p[1][j], universo.p[2][j], universo.m[j], universo, aceleracion, j, Fuerzas, Aceleraciones);
            //cout <<  endl;
            //cout << "Vector aceleracion calculado" << endl;
            //Calculo de las nuevas velocidades y posiciones en los 3 ejes
            for(int k = 0; k < 3; k++){
                newIt.v[k].push_back((universo.v[k][j] + aceleracion[k]*time_step));
                //cout << "Nueva aceleracion " << k << " de asteroide " << j << " " << aceleracion[k] << endl;
                //cout << "Nueva velocidad " << k << " de asteroide " << j << " " << newIt.v[k][j] << endl;
                newIt.p[k].push_back(universo.p[k][j] + newIt.v[k][j]*time_step);
                //cout << "Nueva posicion " << k << " de asteroide " << j << " " << newIt.p[k][j] << endl;
                //cout << "Calculando rebotes" << endl;
                //Calculo de rebotes
                if(newIt.p[k][j] < 0.0){
                    newIt.p[k][j] = 0.0;
                    newIt.v[k][j] *= -1.0;
                }
                if(newIt.p[k][j] > size_enclosure){
                    newIt.p[k][j] = size_enclosure;
                    newIt.v[k][j] *= -1.0;
                }
                //cout << "rebotes calculados" << endl;
            }
            Velocidades << "Velocity(" << j << "): " << newIt.v[0][j] << " " << newIt.v[1][j] << " " << newIt.v[2][j] << "(+ " << newIt.v[0][j] << " " << newIt.v[1][j] << " " << newIt.v[2][j] << ")" << endl;
            //Cargar las masas en la nueva iteracion
            newIt.m.push_back(universo.m[j]);
        }
        //cout << "Actualizando universo" << endl;
        //Actualizacion del universo
        universo = newIt;
        //cout << "Universo actualizado" << endl;
        for (long unsigned int j = 0; j < universo.p[0].size(); j++){
            Posiciones << "Position(" << j << ") = " << universo.p[0][j] << " " << universo.p[1][j] << " " << universo.p[2][j] << endl;
        }
    }

    Fuerzas.close();
    Aceleraciones.close();
    Velocidades.close();
    Posiciones.close();

    //Almacenar configuracion final
    guardarDatos(universo, "final_config.txt", size_enclosure, time_step, int(universo.p[0].size()));
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

    simulation(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atof(argv[4]), atof(argv[5]));

    return 0;
}
