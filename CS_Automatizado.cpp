#include "Murmur/murmurhash32.hpp"
#include <iostream>
#include <vector>
#include <cstdint>
#include <climits>
#include <limits>
#include <seqan/seq_io.h>
using namespace std;

uint64_t canonical_kmer (uint64_t kmer, uint k = 31)
{
    uint64_t reverse = 0;
    uint64_t b_kmer = kmer;

    kmer = ((kmer >> 2)  & 0x3333333333333333UL) | ((kmer & 0x3333333333333333UL) << 2);
    kmer = ((kmer >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((kmer & 0x0F0F0F0F0F0F0F0FUL) << 4);
    kmer = ((kmer >> 8)  & 0x00FF00FF00FF00FFUL) | ((kmer & 0x00FF00FF00FF00FFUL) << 8);
    kmer = ((kmer >> 16) & 0x0000FFFF0000FFFFUL) | ((kmer & 0x0000FFFF0000FFFFUL) << 16);
    kmer = ( kmer >> 32                        ) | ( kmer                         << 32);
    reverse = (((uint64_t)-1) - kmer) >> (8 * sizeof(kmer) - (k << 1));

    return (b_kmer < reverse) ? b_kmer : reverse;
}


class CountSketch {
private:
    int w; // ancho
    int d; // número de hashes
    vector<vector<int64_t>> M;

    int sign_hash(uint64_t elemento, int j) {
        // Hash de signo
        uint32_t h = murmurhash(&elemento, j + 1337);
        return (h & 1) ? 1 : -1;
    }

public:
    CountSketch(int w, int d) : w(w), d(d), M(d, vector<int64_t>(w, 0)) {}

    void insert(uint64_t elemento) {
        for (int j = 0; j < d; j++) {
            uint32_t h_j = murmurhash(&elemento, j) % w;  // hash de índice
            int s_j = sign_hash(elemento, j);             // hash de signo
            M[j][h_j] += s_j;
        }
    }

    int64_t estimate(uint64_t elemento) {
        vector<int64_t> estimates(d);

        for (int j = 0; j < d; j++) {
            uint32_t h_j = murmurhash(&elemento, j) % w;
            int s_j = sign_hash(elemento, j);
            estimates[j] = M[j][h_j] * s_j;
        }

        // Mediana usando sort
        sort(estimates.begin(), estimates.end());
        size_t mid = d / 2;

        if (d % 2 == 1) {
            return estimates[mid];
        } else {
            return (estimates[mid - 1] + estimates[mid]) / 2;
        }
    }
};

// procesar kmers y insertar al cs, similar a como se hace de forma normal

void process_kmers_countsketch(CountSketch &cs, const string &filename, uint k) {
    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, filename.c_str())) {
        cerr << "ERROR: Could not open file " << filename << endl;
        return;
    }

    seqan::CharString id;
    seqan::IupacString seq;

    while (!atEnd(seqFileIn)) {
        try {
            seqan::readRecord(id, seq, seqFileIn);
        } catch (seqan::ParseError &) {
            break;
        }

        uint64_t kmer = 0;
        uint bases = 0;

        for (size_t i = 0; i < length(seq); i++) {
            uint8_t two_bit = 0;

            switch (char(seq[i])) {
                case 'A': case 'a': two_bit = 0; break;
                case 'C': case 'c': two_bit = 1; break;
                case 'G': case 'g': two_bit = 2; break;
                case 'T': case 't': two_bit = 3; break;
                default: // N u otra base
                    kmer = 0; bases = 0;
                    continue;
            }

            kmer = ((kmer << 2) | two_bit) & ((1ULL << (k*2)) - 1);
            bases++;

            if (bases == k) {
                uint64_t canonical = canonical_kmer(kmer, k);
                cs.insert(canonical);
                bases--; // mover ventana
            }
        }
    }

    close(seqFileIn);
}




uint64_t parse_binary_kmer(const std::string &bin_str, uint k) {
    // Convierte la cadena completa de bits en un uint64_t
    uint64_t value = 0;
    for (char c : bin_str) {
        value <<= 1;
        if (c == '1') value |= 1ULL;
    }
    // ignora los 64-2*k bits del inicio
    uint64_t mask = (1ULL << (2 * k)) - 1;
    return value & mask;
}




int main(int argc, char* argv[]) {


    if (argc < 2) {
        std::cerr << "Uso: " << argv[0] << " lista_archivos.txt [-k tamaño_k]\n";
        return 1;
    }

    std::string lista_archivos = argv[1];
    uint k = 31;

    // Leer argumento -k
    int opt;
    while ((opt = getopt(argc, argv, "k:")) != -1) {
        switch (opt) {
            case 'k': k = std::stoi(optarg); break;
            default: break;
        }
    }

    // Abrir archivo con la lista de rutas
    std::ifstream infile(lista_archivos);
    if (!infile.is_open()) {
        std::cerr << "No se pudo abrir el archivo de lista: " << lista_archivos << "\n";
        return 1;
    }

    // parametros
    int w = 2000000;  
    int d = 10;   
    std::string filename;

    //COMPARACION
    while (std::getline(infile, filename)) {
        // limpiar espacios finales
        filename.erase(filename.find_last_not_of(" \r\n\t") + 1);

        if (filename.empty()) continue; // saltar líneas vacías

        cout << "Procesando archivo: " << filename << endl;

        // Procesar kmers
        CountSketch cs(w, d);
        process_kmers_countsketch(cs, filename, k);

        cout << "Inserción de kmers terminada.\n";

        // Buscar la primera aparición de '/'
        size_t pos = filename.find('/');
        if (pos != std::string::npos) {
            filename = filename.substr(pos + 1);
        }

        //std::string ruta   = "Muestra_REALES/K" + std::to_string(k) + "/" + filename + "_" + std::to_string(k) + "_heavy_hitters.txt";
        std::string ruta   = "Total_Reales/K" + std::to_string(k) + "/" + filename + "_" + std::to_string(k) + "_heavy_hitters.txt";
        
        std::string output = "CS_COMPARACION_K_" + std::to_string(k) + "_" + filename + ".txt";

        std::ifstream hhfile(ruta);
        if (!hhfile.is_open()) {
            std::cerr << "No se pudo abrir el archivo de entrada: " << ruta << std::endl;
            continue; // saltar a la siguiente muestra
        }

        std::ofstream outfile(output);
        if (!outfile.is_open()) {
            std::cerr << "No se pudo crear el archivo de salida: " << output << std::endl;
            continue;
        }

        int n;
        hhfile >> n; // la primera línea: cantidad de HH
        outfile << n << "\n";

        std::string bin;
        int freq;
        while (hhfile >> bin >> freq) {
            if (k == 31 && bin.size() > 2)
                bin = bin.substr(2); // quitar "00"
            
            if (k == 21 && bin.size()> 22)
                bin = bin.substr(22); //quitar los 11 "00"

            uint64_t kmer = std::stoull(bin, nullptr, 2);
            int freq_cs = cs.estimate(kmer);

            outfile << bin << " R=" << freq << " CS=" << freq_cs << "\n";
        }

        hhfile.close();
        outfile.close();

        std::cout << "Archivo generado correctamente: " << output << std::endl;
    }

    infile.close();




    //HH CS
    float phi = 0.0000015;

    int contador_hh = 0;

    std::ifstream infile_hh(lista_archivos); 
    if (!infile_hh.is_open()) {
        std::cerr << "No se pudo abrir el archivo de lista.\n";
        return 1;
    }

    while (std::getline(infile_hh, filename)) {
        //limpiar espacios finales
        filename.erase(filename.find_last_not_of(" \r\n\t") + 1);
        if (filename.empty()) continue;

        std::cout << "Procesando archivo HH CS: " << filename << std::endl;

        // Procesar kmers
        CountSketch cs(w, d);
        process_kmers_countsketch(cs, filename, k);
        cout << "Inserción de kmers terminada.\n";

        //Buscar la primera aparición de '/'
        size_t pos = filename.find('/');
        if (pos != std::string::npos) {
            filename = filename.substr(pos + 1);
        }

        std::string ruta_kmers = "Total_Reales/K" + std::to_string(k) + "/" + filename + "_" + std::to_string(k) + "_kmers.txt";
        std::string output_heavy_hitters = "CS_HH_K_" + std::to_string(k) + "_" + filename + ".txt";
        std::string output_kmers = "CS_kmers_K_" + std::to_string(k) + "_" + filename + ".txt";
     
        std::ifstream hhfile(ruta_kmers);
        if (!hhfile.is_open()) {
            std::cerr << "No se pudo abrir el archivo de entrada: " << ruta_kmers << std::endl;
            continue;
        }

        std::ofstream outfile(output_heavy_hitters);
        std::ofstream outfile2(output_kmers);
        if (!outfile.is_open()) {
            std::cerr << "No se pudo crear el archivo de salida: " << output_heavy_hitters << std::endl;
            continue;
        }
        if (!outfile2.is_open()) {
            std::cerr << "No se pudo crear el archivo de salida: " << output_kmers << std::endl;
            continue;
        }

        int n;
        hhfile >> n; // primera línea: cantidad de HH
        outfile2 <<n << "\n";

        float umbral = phi*n;

        std::string bin;
        int freq;
        while (hhfile >> bin >> freq) {
            if (k == 31 && bin.size() > 2)
                bin = bin.substr(2); // quitar "00"
            
            if (k == 21 && bin.size()> 22)
                bin = bin.substr(22); //quitar los 11 "00"

            uint64_t kmer = std::stoull(bin, nullptr, 2);
            int freq_cs = cs.estimate(kmer); // frecuencia de CountSketch
            outfile2 << bin << " R=" << freq << " CS=" << freq_cs << "\n";
            if(freq_cs >= umbral){
                outfile << bin << " R=" << freq << " CS=" << freq_cs << "\n";
                contador_hh++;
            }

        }
        outfile << "N_HH " << contador_hh << "\n";
        
        hhfile.close();
        outfile.close();
        outfile2.close();
        contador_hh = 0;

        std::cout << "Archivo HH CS generado correctamente: " << output_kmers << std::endl;
    }

    infile_hh.close();


    return 0;
}
