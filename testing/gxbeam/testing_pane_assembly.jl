using Plots, StaticArrays, OpenFASTsr, DelimitedFiles, DynamicStallModels, Rotors, LaTeXStrings

DS = DynamicStallModels
of = OpenFASTsr


path = dirname(@__FILE__)
cd(path)



ofpath = "./OpenFAST_NREL5MW_modified" 

inputfile = of.read_inputfile("NREL5MW_input.fst", ofpath)
inflowwind = of.read_inflowwind("NREL5MW_inflowwind.dat", ofpath)
adfile = of.read_adfile("NREL5MW_ADfile.dat", ofpath)
adblade = of.read_adblade("NREL5MW_adblade.dat", ofpath)
edfile = of.read_edfile("NREL5MW_edfile.dat", ofpath)
bdfile = of.read_bdfile("NREL5MW_bdfile.dat", ofpath)
bdblade = of.read_bdblade("NREL5MW_bdblade.dat", ofpath)



assembly = of.make_assembly(edfile, bdfile, bdblade)


outs = Rotors.pane_assembly(assembly; ne=60)



nothing