# Description
This project aims to couple FEBio with an MFEM-based FSI solver [iFSI](https://github.com/ryanblack129/iFSI).

Please see INSTALL.md for instructions on building the plugin. 

## Linking to FEBio
To load the plugin into the FEBio solver, please update your febio.xml file to import the plugin library. For example,

```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_config version="3.0">
	<default_linear_solver type="pardiso"></default_linear_solver>
	<import>/Path/to/plugin/library/libFEBioMFEMiFSI.dylib</import>
</febio_config>
```

## Linking to FEBio Studio
The plugin can be linked to FEBio Studio via the plugin repository. Click the Plugin Repository button in FEBio Studio, then Load Local Plugin and navigate to your plugin library file.

# Notes:
-This repository uses **Git subtree** to include the iFSI code. Please see SUBTREE.md for details.
