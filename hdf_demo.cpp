
#include "hdf_demo.h"

//#define PROCESS_HDF_FILE

int main(int argc, char **argv)
{
    std::string fname = default_filename;
    if(argc > 1)
    {
        fname = argv[1];
    }

    // Uncomment this line to disable the HDF automatic stack printing
    //H5::Exception::dontPrint();

#ifdef PROCESS_HDF_FILE
    std::string processed_filename = default_processed_prefix + fname;
    std::cout << "Processing HDF5 File: " << fname << std::endl;
    std::cout << "Output File: " << processed_filename << std::endl;
    std::cout << std::endl;

    process_hdf_file(fname, processed_filename);
#else
    std::cout << "Reading HDF5 File: " << fname;
    std::cout << std::endl << std::endl;

    // Inspect the contents of the file to discover what data is available
    inspect_hdf_file(fname);

    // Read some of the data in the file
    read_instrument_data(fname);
#endif

    std::cout << std::endl;
    std::cout << "Done" << std::endl;

#ifdef PROCESS_HDF_FILE
    std::cout << "Press Enter to Exit";
    getchar();
#endif
    return 0;
}
