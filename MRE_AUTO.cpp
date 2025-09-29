#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdint>
using namespace std;

// esta cosa procesa un archivo de kmers o hh y devuelve MAE y MRE
bool procesarArchivo(const string &archivo, double &MAE, double &MRE, int &kmerCount) {
    ifstream fin(archivo);
    if (!fin.is_open()) {
        cerr << "No se pudo abrir el archivo: " << archivo << endl;
        return false;
    }

    uint64_t totalKmers;
    fin >> totalKmers;
    string linea;
    getline(fin, linea); // consumir resto de la primera línea

    double sumaErrorRelativo = 0;
    double sumaErrorAbsoluto = 0;
    int contadorKmersUnicos = 0;

    while (getline(fin, linea)) {
        if (linea.empty()) continue;

        stringstream ss(linea);
        string kmer, R_str, TS_str;
        int R, TS;

        ss >> kmer >> R_str >> TS_str;
        R = stoi(R_str.substr(R_str.find('=') + 1));
        TS = stoi(TS_str.substr(TS_str.find('=') + 1));

        sumaErrorAbsoluto += abs(R - TS);
        if (R != 0)
            sumaErrorRelativo += (double)abs(R - TS) / R;
        contadorKmersUnicos++;
    }

    fin.close();

    if (contadorKmersUnicos == 0) return false;

    MAE = sumaErrorAbsoluto / contadorKmersUnicos;
    MRE = sumaErrorRelativo / contadorKmersUnicos;
    kmerCount = contadorKmersUnicos;
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Uso: " << argv[0] << " archivo.txt\n";
        return 1;
    }

    string listaArchivos = argv[1];
    int k = 21; 
    int d = 4; 

    // Parsear argumentos 
    for (int i = 2; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-k" && i + 1 < argc) {
            k = stoi(argv[++i]);
        } else if (arg == "-d" && i + 1 < argc) {
            d = stoi(argv[++i]);
        } else {
            cerr << "Argumento desconocido: " << arg << endl;
            return 1;
        }
    }

    ifstream fin(listaArchivos);
    if (!fin.is_open()) {
        cerr << "No se pudo abrir el archivo de lista: " << listaArchivos << endl;
        return 1;
    }

    string filename;
    double sumaMAE = 0, sumaMRE = 0;
    int archivosProcesados = 0;

    while (getline(fin, filename)) {
        filename.erase(filename.find_last_not_of(" \r\n\t") + 1);
        if (filename.empty()) continue;

        // Quitar prefijo de carpeta si existe
        size_t pos = filename.find('/');
        if (pos != string::npos) {
            filename = filename.substr(pos + 1);
        }

        // Construir ruta final
        //string ruta = "TS_Resultados/Muestra/HH/K" + to_string(k) + "/D=" + to_string(d) + "/TS_kmers_K_" + to_string(k)+  "_" + filename + ".txt";
        string ruta = "TS_Resultados/Total_D=8/HH/K" + to_string(k) + "/FIXED_TS_HH_K_" + to_string(k)+  "_" + filename + ".txt";

        double MAE, MRE;
        int kmerCount;
        if (procesarArchivo(ruta, MAE, MRE, kmerCount)) {
            cout << "Archivo: " << ruta << ", k-mers: " << kmerCount 
                 << ", MAE=" << MAE << ", MRE=" << MRE << endl;
            sumaMAE += MAE;
            sumaMRE += MRE;
            archivosProcesados++;
        } else {
            cerr << "Error procesando archivo: " << ruta << endl;
        }
    }

    fin.close();

    if (archivosProcesados == 0) {
        cerr << "No se procesó ningún archivo correctamente.\n";
        return 1;
    }

    cout << "\nPromedio sobre " << archivosProcesados << " archivos:\n";
    cout << "MAE promedio: " << (sumaMAE / archivosProcesados) << endl;
    cout << "MRE promedio: " << (sumaMRE / archivosProcesados) << endl;

    return 0;
}
