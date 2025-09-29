#include <bitset> // Required for std::bitset
#include <fstream>
#include <omp.h>
#include <iostream>
#include <seqan/seq_io.h>
#include <unordered_set>
#include <unordered_map>

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




void count_kmers(std::unordered_map<uint64_t, uint64_t> &kmer_counts,
                 const std::string &filename, uint k)
{
    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, filename.c_str()))
    {
        std::cerr << "ERROR: Could not open the file " << filename << ".\n";
        return;
    }

    seqan::CharString id;
    seqan::IupacString seq;

    while (!atEnd(seqFileIn))
    {
        try {
            seqan::readRecord(id, seq, seqFileIn);
        }
        catch (seqan::ParseError &)
        {
            break;
        }

        uint64_t kmer = 0;
        uint bases = 0;

        for (size_t i = 0; i < length(seq); ++i)
        //for (size_t i = 0; i < 100; i++)

        {
            uint8_t two_bit = 0;

            switch (char(seq[i]))
            {
                case 'A': case 'a': two_bit = 0; break;
                case 'C': case 'c': two_bit = 1; break;
                case 'G': case 'g': two_bit = 2; break;
                case 'T': case 't': two_bit = 3; break;
                default: // base desconocida (N u otra)
                    kmer = 0;
                    bases = 0;
                    continue; // descarta esta base
            }

            kmer = ((kmer << 2) | two_bit) & ((1ULL << (k*2)) - 1);
            bases++;

            if (bases == k)
            {
                uint64_t canonical = canonical_kmer(kmer, k);
                kmer_counts[canonical]++; // Incrementa la frecuencia
                bases--; // para desplazar ventana de k-mer
            }
        }
    }

    close(seqFileIn);
}


std::unordered_map<uint64_t, uint64_t> heavy_hitters(std::unordered_map<uint64_t, uint64_t> &kmer_counts, float phi, int N){
    std::unordered_map<uint64_t, uint64_t> heavy_hitters_map;
    float threshold = phi * N;

    std::cout << "Threshold (phi * N): " << threshold << "\n";
    for (const auto &pair : kmer_counts) {
        if(pair.second >= threshold) { // si f(e) del kmers, es mayor a phi*N
            heavy_hitters_map[pair.first] = pair.second;
        }
    }
    return heavy_hitters_map;

}




int main(int argc, char *argv[])
{
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






    std::string filename;
    while (std::getline(infile, filename)) {
        filename.erase(filename.find_last_not_of(" \r\n\t") + 1);


        std::cout << "Procesando archivo: " << filename << "\n";


        if (filename.empty()) continue;

        int total_kmers = 0;
        std::unordered_map<uint64_t, uint64_t> kmer_counts;
        count_kmers(kmer_counts, filename, k);

        for (const auto &pair : kmer_counts) {
            total_kmers += pair.second;
        }

        float phi = 0.0000015;
        int N = total_kmers;
        std::unordered_map<uint64_t, uint64_t> heavy_hitters_map =
            heavy_hitters(kmer_counts, phi, N);

        
        std::string base_name = filename.substr(filename.find_last_of("/\\") + 1);

        
        // Guardar los k-mers
        std::string out_kmers = base_name + "_" + std::to_string(k) + "_kmers.txt";
        FILE *output_kmers = fopen(out_kmers.c_str(), "w");
        if (!output_kmers) {
            std::cerr << "No se pudo crear archivo: " << out_kmers << "\n";
            continue;
        }

        fprintf(output_kmers, "%d\n", total_kmers);
        for (const auto &pair : kmer_counts) {
            std::bitset<64> binary_rep(pair.first);
            fprintf(output_kmers, "%s %llu\n",
                    binary_rep.to_string().c_str(),
                    static_cast<unsigned long long>(pair.second));
        }
        fclose(output_kmers);




        // Guardar los heavy hitters
        std::string out_hh = base_name + "_" + std::to_string(k) + "_heavy_hitters.txt";
        FILE *output_hh = fopen(out_hh.c_str(), "w");
        if (!output_hh) {
            std::cerr << "No se pudo crear archivo: " << out_hh << "\n";
            continue;
        }

        fprintf(output_hh, "%zu\n", heavy_hitters_map.size());
        for (const auto &pair : heavy_hitters_map) {
            std::bitset<64> binary_rep(pair.first);
            fprintf(output_hh, "%s %llu\n",
                    binary_rep.to_string().c_str(),
                    static_cast<unsigned long long>(pair.second));
        }
        fclose(output_hh);

        std::cout << "Procesado: " << filename
                  << " -> " << out_kmers
                  << " y " << out_hh << "\n";
    }

    return 0;
}
