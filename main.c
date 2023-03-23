// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

// version 2023/03/23

#include <unistd.h>
#include <groan.h>

void destroy_selections(atom_selection_t **selections, const size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        free(selections[i]);
    }

    free(selections);
}

/*
 * Parses command line arguments.
 * Returns zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments(
        int argc, 
        char **argv,
        char **gro_file,
        char **ndx_file,
        char **output_file,
        char **selection,
        char **phosphate,
        int *empty) 
{
    int gro_specified = 0;

    int opt = 0;
    while((opt = getopt(argc, argv, "c:n:o:s:p:eh")) != -1) {
        switch (opt) {
        // help
        case 'h':
            return 1;
        // gro file to read
        case 'c':
            *gro_file = optarg;
            gro_specified = 1;
            break;
        // ndx file to read
        case 'n':
            *ndx_file = optarg;
            break;
        // output file name
        case 'o':
            *output_file = optarg;
            break;
        // selected atoms
        case 's':
            *selection = optarg;
            break;
        // headgroup identifier
        case 'p':
            *phosphate = optarg;
            break;
        // create empty groups
        case 'e':
            *empty = 1;
            break;
        default:
            //fprintf(stderr, "Unknown command line option: %c.\n", opt);
            return 1;
        }
    }

    if (!gro_specified) {
        fprintf(stderr, "Gro file must always be supplied.\n");
        return 1;
    }
    return 0;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s -c GRO_FILE [OPTION]...\n", program_name);
    printf("\nOPTIONS\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-s STRING        selection of membrane lipids (default: Membrane)\n");
    printf("-p STRING        selection of lipid head identifiers (default: name PO4)\n");
    printf("-o STRING        output ndx file (optional)\n");
    printf("-e               also create empty ndx groups (optional)\n");
    printf("\n");
}

void write_ndx_group(FILE *stream, const char *name, const atom_selection_t *selection)
{
    fprintf(stream, "[ %s ]\n", name);
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        fprintf(stream, "%4ld ", selection->atoms[i]->gmx_atom_number);

        if ((i + 1) % 15 == 0 || i + 1 == selection->n_atoms) fprintf(stream, "\n");
    }
}

/*! @brief Creates ndx groups for lipids distinguishing between membrane leaflets. Returns the number of ndx groups or 0 if no groups were created. */
size_t create_groups(
        const atom_selection_t *membrane, 
        const atom_selection_t *phosphates,
        const list_t *residue_names,
        atom_selection_t ***ndx_groups,
        box_t box) 
{
    // split lipid atoms into individual residues
    atom_selection_t **residues = NULL;
    size_t n_residues = selection_splitbyres(membrane, &residues);

    if (residues == NULL || n_residues == 0) {
        free(residues);
        fprintf(stderr, "Could not split atoms based on residue number.\n");
        return 0;
    }

    // calculate membrane center
    vec_t center = {0.0};
    if (center_of_geometry(membrane, center, box) != 0) {
        fprintf(stderr, "Could not calculate center of geometry for membrane lipids.\n");
        destroy_selections(residues, n_residues);
        return 0;
    }

    // allocate memory for ndx_groups
    size_t n_groups = residue_names->n_items * 2;
    *ndx_groups = calloc(n_groups, sizeof(atom_selection_t *));
    size_t *allocated = calloc(n_groups, sizeof(size_t));

    for (size_t i = 0; i < n_groups; ++i) {
        allocated[i] = 64;
        (*ndx_groups)[i] = selection_create(allocated[i]);
    }

    // loop through all the residues
    for (size_t i = 0; i < n_residues; ++i) {
        char *resname = residues[i]->atoms[0]->residue_name;

        // get lipid phosphate
        atom_selection_t *phosphate = selection_intersect(residues[i], phosphates);
        if (phosphate == NULL || phosphate->n_atoms <= 0) {
            fprintf(stderr, "No phosphate detected for lipid %s (resid %d).\n", resname, residues[i]->atoms[0]->residue_number);
            free(phosphate);
            free(allocated);
            destroy_selections(residues, n_residues);
            destroy_selections(*ndx_groups, n_groups);
            *ndx_groups = NULL;
            return 0;
        }
        if (phosphate->n_atoms > 1) {
            fprintf(stderr, "Multiple phosphates detected for lipid %s (resid %d).\n", resname, residues[i]->atoms[0]->residue_number);
            free(phosphate);
            free(allocated);
            destroy_selections(residues, n_residues);
            destroy_selections(*ndx_groups, n_groups);
            *ndx_groups = NULL;
            return 0;
        }

        // assign lipid into leaflet
        // 1 -> upper, 0 -> lower
        atom_t *pho = phosphate->atoms[0];
        size_t classification = 0;
        if (distance1D(pho->position, center, z, box) > 0) classification = 1;

        // assign lipid into an ndx group
        int index = list_index(residue_names, resname);
        if (index < 0) {
            fprintf(stderr, "Internal Error. Inconsistency in residue names. Residue name %s of resid %d was not found in a list of detected residue names.\n", resname, residues[i]->atoms[0]->residue_number);
            fprintf(stderr, "This should never happen.\n");
            free(phosphate);
            free(allocated);

            destroy_selections(residues, n_residues);
            destroy_selections(*ndx_groups, n_groups);
            *ndx_groups = NULL;
            return 0;
        }

        selection_add(&((*ndx_groups)[2 * index + classification]), &allocated[2 * index + classification], residues[i]);

        free(phosphate);
    }

    free(allocated);
    destroy_selections(residues, n_residues);

    return n_groups;



}

int main(int argc, char **argv)
{
    // get arguments
    char *gro_file = NULL;
    char *ndx_file = "index.ndx";
    char *output_file = NULL;
    char *selected = "Membrane";
    char *phosphate = "name PO4";
    int empty = 0;

    int return_code = 0;

    if (get_arguments(argc, argv, &gro_file, &ndx_file, &output_file, &selected, &phosphate, &empty) != 0) {
        print_usage(argv[0]);
        return 1;
    }
    
    // read gro file
    system_t *system = load_gro(gro_file);
    if (system == NULL) return 1;

    // read ndx file; ignore if this fails
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // select all atoms
    atom_selection_t *all = select_system(system);

    // select membrane lipids
    atom_selection_t *membrane = smart_select(all, selected, ndx_groups);
    if (membrane == NULL) {
        fprintf(stderr, "Could not understand the selection query '%s'.\n", selected);

        dict_destroy(ndx_groups);
        free(all);
        free(system);
        return 1;
    }

    if (membrane->n_atoms == 0) {
        fprintf(stderr, "No membrane lipids ('%s') found.\n", selected);

        dict_destroy(ndx_groups);
        free(membrane);
        free(system);
        free(all);
        return 1;
    }

    // select phosphates
    atom_selection_t *phosphates = smart_select(all, phosphate, ndx_groups);
    if (phosphates == NULL || phosphates->n_atoms == 0) {
        fprintf(stderr, "No phosphates ('%s') found.\n", phosphate);

        dict_destroy(ndx_groups);
        free(phosphates);
        free(membrane);
        free(all);
        free(system);
        return 1;
    }

    // get residue names
    list_t *residue_names = selection_getresnames(membrane);

    // create new ndx groups
    atom_selection_t **lipids_leaflets = NULL;
    size_t n_groups = 0;
    if ( (n_groups = create_groups(membrane, phosphates, residue_names, &lipids_leaflets, system->box)) == 0) {
        fprintf(stderr, "Failed to create ndx groups.\n");
        return_code = 1;
        goto main_end;
    }

    // open the output file
    FILE *output = NULL;
    FILE *test = NULL;
    if (output_file == NULL) output = stdout;
    // if the file exists, append
    else if ((test = fopen(output_file, "r")) != NULL) {
        fclose(test);
        output = fopen(output_file, "a");
    } else {
        output = fopen(output_file, "w");
    }

    if (output == NULL) {
        fprintf(stderr, "The output ndx file could not be opened.\n");
        return_code = 1;
        goto main_end;
    }

    // write out the lipids_leaflets ndx groups
    for (size_t i = 0; i < n_groups; ++i) {

        if (!empty && lipids_leaflets[i]->n_atoms == 0) continue;

        char group_name[100] = "";
        char *resname = list_get(residue_names, i / 2);
        if (resname == NULL) {
            fprintf(stderr, "Internal error. Reaching element of index %ld in a list_t of length %ld", i / 2, residue_names->n_items);
            fprintf(stderr, "This should never happen.\n");
            return_code = 1;
            goto main_end;
        }

        strncpy(group_name, resname, 99);
        if (i % 2 == 0) {
            strcat(group_name, "_lower");
        } else {
            strcat(group_name, "_upper");
        }

        write_ndx_group(output, group_name, lipids_leaflets[i]);

    }

    if (output != stdout) fclose(output);

    main_end:
    list_destroy(residue_names);
    destroy_selections(lipids_leaflets, n_groups);
    dict_destroy(ndx_groups);
    free(phosphates);
    free(membrane);
    free(all);
    free(system);

    return return_code;
}