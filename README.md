# Description
This project aims to couple FEBio with an MFEM-based FSI solver ([iFSI](https://github.com/ryanblack129/iFSI)).

Please see INSTALL.md for instructions on building the plugin. 

## Linking to FEBio
To load the plugin into the FEBio solver, please update your febio.xml file to import the plugin library.

## Linking to FEBio Studio
The plugin can then be linked to FEBio Studio via the plugin repository. Click the Plugin Repository button in FEBio Studio, then Load Local Plugin and point to the plugin library.

# Notes:
-This repository uses **Git subtree** to include the iFSI code. Please see SUBTREE.md for details.
