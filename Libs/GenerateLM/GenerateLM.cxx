/*=========================================================================

  Portions (c) Copyright 2006 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   GenerateLM
  Module:    $URL: http://svn.na-mic.org:8000/svn/NAMICSandBox/trunk/CommandLineAPI/GenerateCLP.cxx $
  Date:      $Date: 2006-04-21 15:51:08 -0400 (Fri, 21 Apr 2006) $
  Version:   $Revision: 957 $

=========================================================================*/

/* Generate loadable mode support from an xml description
   Usage: GenerateLM input_xml_file output_include_file ouput_source_file

   This program generates source code that allows module detection at
   run-time.  The module description is in an xml file.  The output
   files defines the required entry points that will be called during
   startup to initialize the module.

   Typical usage is:
   GenerateLM foo.xml fooLib.h fooLib.cxx

   fooLib.cxx contains:
   #include fooLib.h
*/


#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include "expat.h"
#include <string>
#include <vector>
#include <itksys/SystemTools.hxx>

// use GenerateCLP to make the header file that implements PARSE_ARGS
#include "GenerateLMCLP.h"

#include "LoadableModuleDescriptionParser.h"
#include "LoadableModuleDescription.h"

/* Generate the preamble to the code. This includes the required
 * include files and code to process comma separated arguments.
 */
void GeneratePre(std::ofstream &, LoadableModuleDescription &, int, char *[]);

/* Generate the code that echos the XML file that describes the
 * command line arguments.
 */
void GenerateSourceIncludes(std::ofstream &, const std::string&, int, char *[]);
void GenerateExports(std::ofstream &);
void GenerateModuleDataSymbols(std::ofstream &, const std::string&, const std::string&);
void GenerateModuleEntryPoints(std::ofstream &);

int
main(int argc, char *argv[])
{
  PARSE_ARGS;
  LoadableModuleDescription module;
  LoadableModuleDescriptionParser parser;

  // Read the XML file
  std::ifstream fin(InputXML.c_str(),std::ios::in|std::ios::binary);
  if (fin.fail())
    {
    std::cerr << argv[0] << ": Cannot open " << InputXML << " for input" << std::endl;
    perror(argv[0]);
    return EXIT_FAILURE;
    }

  // Get the length of the file
  fin.seekg (0, std::ios::end);
  size_t len = fin.tellg();
  fin.seekg (0, std::ios::beg);
  char * XML = new char[len+1];
  fin.read (XML, len);
  XML[len] = '\0';

  // Parse the module description
  std::cerr << "GenerateLM";
  for (int i = 1; i < argc; i++)
    {
    std::cerr << " " << argv[i];
    }
  std::cerr << std::endl;

  if (parser.Parse(XML, module))
    {
    std::cerr << "GenerateLM: One or more errors detected. Code generation aborted." << std::endl;
    return EXIT_FAILURE;
    }

  // Do the hard stuff

  {
    std::ofstream sout(OutputCxx.c_str(),std::ios::out);
    if (sout.fail())
      {
        std::cerr << argv[0] << ": Cannot open " << OutputCxx << " for output" << std::endl;
        perror(argv[0]);
        return EXIT_FAILURE;
      }
    
    // make source cxx
    std::string fname = itksys::SystemTools::GetFilenameName(InputXML);

    std::string classname = fname.substr(0, fname.rfind('.'));

    GenerateSourceIncludes(sout, classname, argc, argv);
    GenerateModuleDataSymbols(sout, InputXML, classname);
    sout.close();
  }

  {
    std::ofstream sout(OutputHeader.c_str(),std::ios::out);
    if (sout.fail())
      {
        std::cerr << argv[0] << ": Cannot open " << OutputHeader << " for output" << std::endl;
        perror(argv[0]);
        return EXIT_FAILURE;
      }


    // make header


    GeneratePre(sout, module, argc, argv);
    GenerateExports(sout);
    GenerateModuleEntryPoints(sout);
    sout.close();
  }

  return (EXIT_SUCCESS);
}

void GeneratePre(std::ofstream &sout, LoadableModuleDescription &module, int argc, char *argv[])
{
  sout << "// This file was automatically generated by:" << std::endl;
  sout << "// ";
  for (int i = 0; i < argc; i++)
    {
    sout << " " << argv[i];
    }
  sout << std::endl;
  sout << "//" << std::endl;
  sout << "#include <stdio.h>" << std::endl;
  sout << "#include <stdlib.h>" << std::endl;
  sout << "#include <iostream>" << std::endl;
  sout << "#include <string.h>" << std::endl;
  sout << std::endl;
}

void GenerateSourceIncludes(std::ofstream &sout, const std::string &classname, int argc, char *argv[])
{
  sout << "// This file was automatically generated by:" << std::endl;
  sout << "// ";
  for (int i = 0; i < argc; i++)
    {
    sout << " " << argv[i];
    }
  sout << std::endl;
  sout << "#include \"" << classname << "Lib.h\"" << std::endl;
  sout << std::endl;
  sout << "#include \"vtk" << classname << "Logic.h\"" << std::endl;
  sout << "#include \"vtk" << classname << "GUI.h\"" << std::endl;
  sout << std::endl;
}

void GenerateExports(std::ofstream &sout)
{
  sout << "#ifdef WIN32" << std::endl;
  sout << "#define Module_EXPORT __declspec(dllexport)" << std::endl;
  sout << "#else" << std::endl;
  sout << "#define Module_EXPORT " << std::endl;
  sout << "#endif" << std::endl;
  sout << std::endl;
}

void GenerateModuleDataSymbols(std::ofstream &sout, const std::string &XMLFile, const std::string &classname)
{
  sout << "extern \"C\" {" << std::endl;
  sout << "Module_EXPORT char " << classname << "ModuleDescription[] = " << std::endl;

  std::string line;
  std::ifstream fin(XMLFile.c_str(),std::ios::in);
  while (!fin.eof())
    {
    std::getline( fin, line );

    // replace quotes with escaped quotes
    std::string cleanLine;
    for (size_t j = 0; j < line.length(); j++)
      {
      if (line[j] == '\"')
        {
        cleanLine.append("\\\"");
        }
      else
        {
        cleanLine.append(1,line[j]);
        }
      }
    sout << "\"" << cleanLine << "\\n\"" << std::endl;
    }
  sout << ";" << std::endl;
  sout << "}" << std::endl;
  sout << std::endl;

  fin.close();

  sout << "char*" << std::endl;
  sout << "GetLoadableModuleDescription()" << std::endl;
  sout << "{" << std::endl;
  sout << "  return " << classname << "ModuleDescription;" << std::endl;
  sout << "}" << std::endl;
  sout << std::endl;
  sout << "void*" << std::endl;
  sout << "GetLoadableModuleGUI()" << std::endl;
  sout << "{" << std::endl;
  sout << "  return vtk" << classname << "GUI::New ( );" << std::endl;
  sout << "}" << std::endl;
  sout << std::endl;
  sout << std::endl;
  sout << "void*" << std::endl;
  sout << "GetLoadableModuleLogic()" << std::endl;
  sout << "{" << std::endl;
  sout << "  return vtk" << classname << "Logic::New ( );" << std::endl;
  sout << "}" << std::endl;

}

void GenerateModuleEntryPoints(std::ofstream &sout)
{
  sout << "extern \"C\" {" << std::endl;
  sout << "  Module_EXPORT char* GetLoadableModuleDescription();" << std::endl;
  sout << "  Module_EXPORT void* GetLoadableModuleGUI();" << std::endl;
  sout << "  Module_EXPORT void* GetLoadableModuleLogic();" << std::endl;
  sout << "}" << std::endl;
  sout << std::endl;
}
