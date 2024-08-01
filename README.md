---
title: "README SpeSpeNet"
author: "Abraham Lucas van Eijnatten"
date: "2024-07-17"
---

## Abstract

Correlation networks are commonly used to explore microbiome data. In these networks, nodes are taxa and edges represent correlations between their abundance patterns across samples. As clusters of correlating taxa (co-abundance clusters) often  indicate a shared response to environmental drivers, network visualization contributes to system understanding. Currently, most tools for creating and visualizing co-abundance networks from microbiome data either require the researcher to have coding skills, or they are not user-friendly, with high time expenditure and limited customizability. Furthermore, existing tools lack focus on the relationship between environmental drivers and the structure of the microbiome, even though many edges in correlation networks can be understood through a shared relationship of two taxa with the environment. For these reasons we developed SpeSpeNet (Species-Species Network, https://tbb.bio.uu.nl/SpeSpeNet), a practical and user-friendly R-shiny tool to construct and visualize correlation networks from taxonomic abundance  tables. The details of data preprocessing, network construction, and visualization are automated, require no programming ability for the web version, and are highly customizable, including associations with user-provided environmental data.

Find the SpeSpeNet tool here: https://tbb.bio.uu.nl/SpeSpeNet. 

Find the paper presenting SpeSpeNet here:

## Description files on this github:

Source code: map containing three scripts. App.R contains the code for the ui and server of SpeSpeNet. nw_functions.R contains functions sourced by app.R. Launch_tool.R is the script that loads in the dependencies and launches the desktop version of the tool. 

Example_data: Map containing three data sets that are in SpeSpeNet compatible format and can be used to try SpeSpeNet. The three datasets correspond to the three case studies presented in the paper (Hauptfeld et al, 2022; Deutschmann et al, 2023; Brenzinger et al, 2021).

Figure_scripts: Map with the scripts to make the figures of the SpeSpeNet paper. These scripts are a linearized version of the SpeSpeNet source code.

SpeSpeNet_manual.docx: In depth description of how to use SpeSpeNet.  

## Remarks

Currently it is very difficult to obtain a workable installation of R to run SpeSpeNet due to the MGnifyR package no longer being compatible with the other packages. Hence the webtool is recommended.

## Contact

For issues, requests or discussions please open an issue on the github or contact [a.l.vaneijnatten\@uu.nl](mailto:a.l.vaneijnatten@uu.nl) or [l.b.snoek\@uu.nl](mailto:l.b.snoek@uu.nl)
