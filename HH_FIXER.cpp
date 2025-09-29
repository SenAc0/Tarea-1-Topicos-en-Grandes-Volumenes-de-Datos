#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

// procesar archivo individual
void procesarArchivo(const string &archivo) {
    ifstream fin(archivo);
    if (!fin.is_open()) {
        cerr << "No se pudo abrir: " << archivo << "\n";
        return;
    }

    vector<string> lineas;
    string linea;

    while (getline(fin, linea)) {
        // Limpiar posibles caracteres sucios al final
        linea.erase(linea.find_last_not_of(" \r\n\t") + 1);
        if (!linea.empty()) {
            lineas.push_back(linea);
        }
    }
    fin.close();

    if (lineas.empty()) return;

    
    string ultima = lineas.back();
    lineas.pop_back();
    if (ultima.rfind("N_HH", 0) == 0) {
        ultima = ultima.substr(5); // quitar "N_HH " (5 caracteres)
        ultima.erase(0, ultima.find_first_not_of(" \t")); // quitar espacios al inicio
    }

    // Insertar al inicio
    lineas.insert(lineas.begin(), ultima);

    // Crear nuevo archivo con prefijo FIXED_
    string nuevoArchivo = "FIXED_" + archivo;
    ofstream fout(nuevoArchivo);
    if (!fout.is_open()) {
        cerr << "No se pudo crear: " << nuevoArchivo << "\n";
        return;
    }

    for (const auto &l : lineas) {
        fout << l << "\n";
    }
    fout.close();

    cout << "Archivo procesado: " << nuevoArchivo << "\n";
}

int main() {
    // Abrir lista de archivos
    ifstream lista("lista_archivos.txt");
    if (!lista.is_open()) {
        cerr << "No se pudo abrir lista_archivos.txt\n";
        return 1;
    }

    string filename;
    while (getline(lista, filename)) {
        // Limpiar caracteres sucios al final y espacios
        filename.erase(filename.find_last_not_of(" \r\n\t") + 1);
        if (filename.empty()) continue;

        procesarArchivo(filename);
    }

    return 0;
}
