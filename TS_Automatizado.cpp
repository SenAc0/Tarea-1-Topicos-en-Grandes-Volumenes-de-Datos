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

class TowerSketch {
private:
    int levels;
    vector<int> widths;
    vector<int> bits;
    vector<uint64_t> maxVals;
    int d;

    // Cada nivel puede estar en uint8_t, uint16_t, uint32_t o uint64_t
    vector<vector<vector<uint8_t>>>  towers8;
    vector<vector<vector<uint16_t>>> towers16;
    vector<vector<vector<uint32_t>>> towers32;
    vector<vector<vector<uint64_t>>> towers64;

public:
    TowerSketch(const vector<int> &widths, const vector<int> &bits, int d)
        : widths(widths), bits(bits), d(d) {
        levels = widths.size();
        maxVals.resize(levels);

        for (int l = 0; l < levels; l++) {
            if (bits[l] <= 8) {
                towers8.resize(levels);
                towers8[l].resize(d, vector<uint8_t>(widths[l], 0));
                maxVals[l] = (1u << bits[l]) - 1;
            } else if (bits[l] <= 16) {
                towers16.resize(levels);
                towers16[l].resize(d, vector<uint16_t>(widths[l], 0));
                maxVals[l] = (1u << bits[l]) - 1;
            } else if (bits[l] <= 32) {
                towers32.resize(levels);
                towers32[l].resize(d, vector<uint32_t>(widths[l], 0));
                maxVals[l] = (bits[l] >= 32) ? UINT32_MAX : ((1u << bits[l]) - 1);
            } else { // 64 bits
                towers64.resize(levels);
                towers64[l].resize(d, vector<uint64_t>(widths[l], 0));
                maxVals[l] = (bits[l] >= 64) ? UINT64_MAX : ((1ULL << bits[l]) - 1);
            }
        }
    }

    void insert(uint64_t elemento) {
        for (int l = 0; l < levels; l++) {
            uint64_t min_val = UINT64_MAX;
            vector<int> h_idx(d);

            // calcular hashes
            for (int j = 0; j < d; j++) {
                h_idx[j] = murmurhash(&elemento, j) % widths[l];
                uint64_t val;
                if (bits[l] <= 8) val = towers8[l][j][h_idx[j]];
                else if (bits[l] <= 16) val = towers16[l][j][h_idx[j]];
                else if (bits[l] <= 32) val = towers32[l][j][h_idx[j]];
                else val = towers64[l][j][h_idx[j]];
                if (val < min_val) min_val = val;
            }

            // conservative Update
            for (int j = 0; j < d; j++) {
                if (bits[l] <= 8) {
                    auto &bucket = towers8[l][j][h_idx[j]];
                    if (bucket == min_val && bucket < maxVals[l]) bucket++;
                } else if (bits[l] <= 16) {
                    auto &bucket = towers16[l][j][h_idx[j]];
                    if (bucket == min_val && bucket < maxVals[l]) bucket++;
                } else if (bits[l] <= 32) {
                    auto &bucket = towers32[l][j][h_idx[j]];
                    if (bucket == min_val && bucket < maxVals[l]) bucket++;
                } else {
                    auto &bucket = towers64[l][j][h_idx[j]];
                    if (bucket == min_val && bucket < maxVals[l]) bucket++;
                }
            }
        }
    }

    uint64_t frecuencia(uint64_t elemento) {
        for (int l = 0; l < levels; l++) {
            uint64_t min_row = UINT64_MAX;
            bool full = true; //  lleno

            for (int j = 0; j < d; j++) {
                int h = murmurhash(&elemento, j) % widths[l];
                uint64_t val;
                if (bits[l] <= 8) val = towers8[l][j][h];
                else if (bits[l] <= 16) val = towers16[l][j][h];
                else if (bits[l] <= 32) val = towers32[l][j][h];
                else val = towers64[l][j][h];

                // Si alguno no está lleno, este nivel sirve
                if (val < maxVals[l]) full = false;

                if (val < min_row) min_row = val;
            }

            // si este nivel no está lleno → devolver este valor
            if (!full) return min_row;
        }

        // si todos los niveles están llenos → devolver el último
        return maxVals.back();
    }
};


// procesar kmers y insertar al ts, similar a como se hace de forma normal
void process_kmers_towersketch(TowerSketch &ts, const string &filename, uint k) {
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
                ts.insert(canonical);
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
    uint k = 21;

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




    

    //configracion 2
    // vector<int> widths = {1500000, 750000, 375000, 150000};
    // vector<int> bits   = {8, 16, 32, 64};
    // int d = 6;

    
    //configuracion 3
    vector<int> widths = {1800000, 900000, 450000, 200000};
    vector<int> bits   = {8, 16, 32, 64};
    int d = 8;
    std::string filename;

    //COMPARACION
    while (std::getline(infile, filename)) {
        // limpiar espacios finales
        filename.erase(filename.find_last_not_of(" \r\n\t") + 1);

        if (filename.empty()) continue; // saltar líneas vacías

        cout << "Procesando archivo: " << filename << endl;

        // Procesar kmers
        TowerSketch ts(widths, bits, d);
        process_kmers_towersketch(ts, filename, k);

        cout << "Inserción de kmers terminada.\n";

        // Buscar la primera aparición de '/'
        size_t pos = filename.find('/');
        if (pos != std::string::npos) {
            filename = filename.substr(pos + 1);
        }

        std::string ruta   = "Muestra_REALES/K" + std::to_string(k) + "/" + filename + "_" + std::to_string(k) + "_heavy_hitters.txt";
        std::string output = "TS_COMPARACION_K_" + std::to_string(k) + "_" + filename + ".txt";

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
            
            if (k==21 && bin.size()> 22)
                bin = bin.substr(22); //quitar los 11 "00"

            uint64_t kmer = std::stoull(bin, nullptr, 2);
            int freq_ts = ts.frecuencia(kmer);

            outfile << bin << " R=" << freq << " TS=" << freq_ts << "\n";
        }

        hhfile.close();
        outfile.close();

        std::cout << "Archivo generado correctamente: " << output << std::endl;
    }

    infile.close();




    //HH TS
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

        std::cout << "Procesando archivo HH TS: " << filename << std::endl;

        // Procesar kmers
        TowerSketch ts(widths, bits, d);
        process_kmers_towersketch(ts, filename, k);
        cout << "Inserción de kmers terminada.\n";

        //Buscar la primera aparición de '/'
        size_t pos = filename.find('/');
        if (pos != std::string::npos) {
            filename = filename.substr(pos + 1);
        }

        std::string ruta_kmers = "Muestra_REALES/K" + std::to_string(k) + "/" + filename + "_" + std::to_string(k) + "_kmers.txt";
        std::string output_heavy_hitters = "TS_HH_K_" + std::to_string(k) + "_" + filename + ".txt";
        std::string output_kmers = "TS_kmers_K_" + std::to_string(k) + "_" + filename + ".txt";
     
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
            
            if (k==21 && bin.size()> 22)
                bin = bin.substr(22); //quitar los 11 "00"

            uint64_t kmer = std::stoull(bin, nullptr, 2);
            int freq_ts = ts.frecuencia(kmer); // frecuencia de TowerSketch
            outfile2 << bin << " R=" << freq << " TS=" << freq_ts << "\n";
            if(freq_ts >= umbral){
                outfile << bin << " R=" << freq << " TS=" << freq_ts << "\n";
                contador_hh++;
            }

        }
        outfile << "N_HH " << contador_hh << "\n";
        
        hhfile.close();
        outfile.close();
        outfile2.close();
        contador_hh = 0;

        std::cout << "Archivo HH TS generado correctamente: " << output_kmers << std::endl;
    }

    infile_hh.close();


    return 0;
}
