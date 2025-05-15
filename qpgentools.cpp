#include "qpgentools.h"
#include "qgenlib/commands.h"
#include "qgenlib/qgen_utils.h"

int32_t test_qpgen(int32_t argc, char** argv);
int32_t match_plp_geno(int32_t argc, char** argv);
int32_t join_plp_geno(int32_t argc, char** argv);
int32_t cmd_geno2tsv(int32_t argc, char** argv);
 
int32_t main(int32_t argc, char** argv) {
  commandList cl;

  BEGIN_LONG_COMMANDS(longCommandlines)
    LONG_COMMAND_GROUP("Utilities for multi-omics", NULL)
    LONG_COMMAND("match-plp-geno", &match_plp_geno, "Check the concordance between pileup and genotypes")
    LONG_COMMAND("join-plp-geno",  &join_plp_geno, "Join the pileup and individual genotypes")

    LONG_COMMAND_GROUP("Utilities for genotype extraction", NULL)
    LONG_COMMAND("geno2tsv", &cmd_geno2tsv, "Extract genotypes into TSV format")

    LONG_COMMAND_GROUP("Other Utilities", NULL)
    LONG_COMMAND("test-qpgen", &test_qpgen, "Test software to check the qpgen library")
  END_LONG_COMMANDS();

  cl.Add(new longCommands("Available Commands", longCommandlines));
  
  if ( argc < 2 ) {
    printf("[qpgentools] -- Another command line toolkit for PLINK files\n\n");
    fprintf(stderr, " Copyright (c) 2023-2025 by Hyun Min Kang\n");
    fprintf(stderr, " Licensed under the Apache License v2.0 and LGPLv3 (see LICENSE for details)\n");    
    fprintf(stderr, "To run a specific command      : %s [command] [options]\n",argv[0]);
    fprintf(stderr, "For detailed instructions, run : %s --help\n",argv[0]);        
    cl.Status();
    return 1;
  }
  else {
    if ( strcmp(argv[1],"--help") == 0 ) {
      cl.HelpMessage();
    }
    else {
      return cl.Read(argc, argv);
    }
  }
  return 0;
}
