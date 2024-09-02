#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Removes water molecules ("HOH") from a PDB file.
 *
 * @param input_file_path Path to the PDB file to be processed.
 */
void remove_water_from_file(const char *input_file_path) {
    char temp_file_path[256];
    snprintf(temp_file_path, sizeof(temp_file_path), "%s.tmp", input_file_path);

    FILE *f_in = fopen(input_file_path, "r");
    FILE *f_out = fopen(temp_file_path, "w");

    if (f_in == NULL || f_out == NULL) {
        fprintf(stderr, "Failed to open file: %s\n", input_file_path);
        if (f_in) fclose(f_in);
        if (f_out) fclose(f_out);
        return;
    }

    char line[256];
    int has_water = 0;
    while (fgets(line, sizeof(line), f_in)) {
        if (strncmp(line, "HETATM", 6) == 0 && strstr(line, "HOH") != NULL) {
            has_water = 1;
        } else {
            fputs(line, f_out);
        }
    }

    fclose(f_in);
    fclose(f_out);

    if (has_water) {
        if (remove(input_file_path) != 0 || rename(temp_file_path, input_file_path) != 0) {
            fprintf(stderr, "Error replacing original file: %s\n", input_file_path);
        }
    } else {
        remove(temp_file_path);  // Remove temporary file if no changes were made
    }
}

/**
 * Main function to handle command line arguments and process PDB files.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 * @return Exit status.
 */
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <pdb_file1> [<pdb_file2> ...]\n", argv[0]);
        return 1;
    }

    for (int i = 1; i < argc; i++) {
        remove_water_from_file(argv[i]);
    }

    return 0;
}
