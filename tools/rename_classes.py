#!/usr/bin/env python3

from fileinput import FileInput
import os

# Repos to process
repos = ["saber", "soca"]

# Replacements
repl = {}
repl["SaberBlockChainBase"] = "BlockChainBase"
repl["SaberBlockChainMaker"] = "BlockChainMaker"
repl["SaberBlockChainFactory"] = "BlockChainFactory"
repl["SaberBlockParametersBase"] = "BlockParametersBase"
repl["SaberCentralBlockBase"] = "CentralBlockBase"
repl["SaberCentralBlockMaker"] = "CentralBlockMaker"
repl["SaberCentralBlockFactory"] = "CentralBlockFactory"
repl["SaberCentralBlockParametersWrapper"] = "CentralBlockParametersWrapper"
repl["SaberCentralBlock"] = "CentralBlockWrapper"
repl["SaberEnsembleBlockChain"] = "EnsembleBlockChain"
repl["SaberGSIBlockChain"] = "GSIBlockChain"
repl["SaberHybridBlockChain"] = "HybridBlockChain"
repl["SaberOuterBlockBase"] = "OuterBlockBase"
repl["SaberOuterBlockChain"] = "OuterBlockChain"
repl["SaberOuterBlockMaker"] = "OuterBlockMaker"
repl["SaberOuterBlockFactory"] = "OuterBlockFactory"
repl["SaberOuterBlockParametersWrapper"] = "OuterBlockParametersWrapper"
repl["SaberOuterBlock"] = "OuterBlock"
repl["SaberParametricBlockChain"] = "ParametricBlockChain"
repl["innerSaberOuterBlockParams"] = "innerOuterBlockParams"
repl["makerSaberDiffusion_"] = "makerDiffusion_"
repl["makerSaberDiffusionFilter_"] = "makerDiffusionFilter_"

# Patterns to check for indentation
patterns = ["const", "atlas::", "eckit::", "oops::", "bool", "int", "float", "double", "std::string","std::shared_ptr"]

# Function to process a single file
def process(filePath):
  # Replacements
  with FileInput(filePath, inplace=True) as file:
    checkNextLine = False
    for line in file:
      if checkNextLine:
        foundPattern = False
        for pattern in patterns:
          if "     " + pattern in line:
            foundPattern = True
            line = line.replace("     " + pattern, pattern)
        if not foundPattern:
          checkNextLine = False
      for key in repl.keys():
        if key in line:
          if not checkNextLine:
            checkNextLine = (key + "(" in line) or (key + "> create(" in line) or (key + "> make(" in line)
          line = line.replace(key, repl[key])
      print(line, end="")

  # Renaming
  for key in repl.keys():
    if key in filePath:
      newFilePath = filePath.replace(key, repl[key])
      os.rename(filePath, newFilePath)
      filePath = newFilePath

# Get SABER tools directory
scriptDir = os.path.dirname(os.path.realpath(__file__))

# Loop over SABER and SOCA directories
for repo in repos:
  # Get source directory
  srcDir = scriptDir + "/../../" + repo + "/src/" + repo

  # Walk through files
  for root, subDirs, fileNames in os.walk(srcDir):
    for fileName in fileNames:
      print("Processing: " + root + "/" + fileName)
      process(root + "/" + fileName)
