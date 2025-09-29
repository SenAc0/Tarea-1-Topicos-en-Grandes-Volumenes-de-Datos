#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <string>
using namespace std;

std::string cleanKmer(const std::string &kmer, int k) {
    std::string bin = kmer;
    if (k == 31 && bin.size() > 2) {
        bin = bin.substr(2); // quitar los 2 ceros sobrantes
    } else if (k == 21 && bin.size() > 22) {
        bin = bin.substr(22); // quitar los 22 ceros sobrantes
    }
    return bin;
}


// heavy hitters reales
bool leerReales(const string &archivo, unordered_set<string> &hhReales, int k) {
    ifstream fin(archivo);
    if (!fin.is_open()) { cerr << "No se pudo abrir archivo real: " << archivo << endl; return false; }

    uint64_t total;
    fin >> total;
    string linea;
    getline(fin, linea); // consumir resto de primera línea

    while (getline(fin, linea)) {
        linea.erase(linea.find_last_not_of(" \r\n\t") + 1);
        if (linea.empty()) continue;
        stringstream ss(linea);
        string kmer;
        int freq;
        if (!(ss >> kmer >> freq)) continue;
        hhReales.insert(cleanKmer(kmer,k));
    }
    fin.close();
    return true;
}

// heavy hitters estimados
bool leerEstimados(const string &archivo, unordered_set<string> &hhEstimados, int k) {
    ifstream fin(archivo);
    if (!fin.is_open()) { cerr << "No se pudo abrir archivo estimado: " << archivo << endl; return false; }

    string linea;
    while (getline(fin, linea)) {
        linea.erase(linea.find_last_not_of(" \r\n\t") + 1);
        if (linea.empty() || linea.rfind("N_HH", 0) == 0) continue;
        stringstream ss(linea);
        string kmer, R_str, TS_str;
        if (!(ss >> kmer >> R_str >> TS_str)) continue;
        hhEstimados.insert(kmer);
    }
    fin.close();
    return true;
}

// Calcular Precision, Recall y F1
void calcularPRF(const unordered_set<string> &reales, const unordered_set<string> &estimados,
                 double &precision, double &recall, double &f1) {
    int TP = 0, FP = 0, FN = 0;

    for (const auto &kmer : estimados)
        if (reales.count(kmer)) TP++;
        else FP++;

    for (const auto &kmer : reales)
        if (!estimados.count(kmer)) FN++;

    precision = (TP + FP > 0) ? double(TP) / (TP + FP) : 0;
    recall = (TP + FN > 0) ? double(TP) / (TP + FN) : 0;
    f1 = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0;
}




int main(int argc, char* argv[]) {
    if (argc < 2) { cerr << "Uso: " << argv[0] << " genomas_muestra.txt [-k 21] [-d 3]\n"; return 1; }

    string listaArchivos = argv[1];
    int k = 21, d = 3;

    for (int i = 2; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-k" && i + 1 < argc) k = stoi(argv[++i]);
        else if (arg == "-d" && i + 1 < argc) d = stoi(argv[++i]);
    }

    ifstream fin(listaArchivos);
    if (!fin.is_open()) { cerr << "No se pudo abrir " << listaArchivos << endl; return 1; }

    string filename;
    double sumaPrecision=0, sumaRecall=0, sumaF1=0;
    int archivosProcesados=0;

    while (getline(fin, filename)) {
        filename.erase(filename.find_last_not_of(" \r\n\t")+1);
        if (filename.empty()) continue;

        // Quitar cualquier prefijo de carpeta para tomar solo el nombre del archivo (fue por lso txts)
        size_t pos = filename.find_last_of('/');
        string baseNombre = (pos != string::npos) ? filename.substr(pos + 1) : filename;

        string archivoReales = "Muestra_REALES/K" + to_string(k) + "/" + baseNombre + "_" + to_string(k) + "_heavy_hitters.txt";
        string archivoEstimados = "resultados TS/muestreo/k" + to_string(k) + "/configuracion " + to_string(d) +
                                "/TS_HH_K_" + to_string(k) + "_" + baseNombre + ".txt";
        
        
        unordered_set<string> hhReales, hhEstimados;

        
        if (!leerReales(archivoReales, hhReales, k) || !leerEstimados(archivoEstimados, hhEstimados, k))
            continue;

        double precision, recall, f1;
        calcularPRF(hhReales, hhEstimados, precision, recall, f1);

         cout << "Archivo reales: " << archivoReales << ", HH leídos: " << hhReales.size() << endl;
        cout << "Archivo estimados: " << archivoEstimados << ", HH leídos: " << hhEstimados.size() << endl;

    
        cout << "Archivo: " << filename
             << ", Precision=" << precision
             << ", Recall=" << recall
             << ", F1=" << f1 << endl;

        sumaPrecision += precision;
        sumaRecall += recall;
        sumaF1 += f1;
        archivosProcesados++;
    }
   

    if (archivosProcesados>0) {
        cout << "\nPromedio sobre " << archivosProcesados << " archivos:\n";
        cout << "Precision promedio: " << sumaPrecision / archivosProcesados << endl;
        cout << "Recall promedio: " << sumaRecall / archivosProcesados << endl;
        cout << "F1-score promedio: " << sumaF1 / archivosProcesados << endl;
    } else {
        cerr << "No se procesó ningún archivo.\n";
    }

    return 0;
}
