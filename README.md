---
title: "README SPESPENET"
author: "Abraham Lucas van Eijnatten"
date: "2024-05-29"
output: html_document
---

## Abstract

Correlation networks are commonly used to explore microbiome data. In these networks, nodes are taxa and edges represent correlations between their relative abundance patterns across samples. Clusters of correlating taxa (co-response clusters) are often indicative of a shared response to environmental drivers. Currently, most tools for network visualizations of microbiome data either require the researcher to have coding skills or are not user-friendly, with high time expenditure and limited customizability. Furthermore, existing tools lack focus on the relationship between environmental drivers and the structure of the microbiome, even though many edges in correlation networks can be understood through a shared relationship of two taxa with the environment. For these reasons we developed SpeSpeNet (Species-Species Network), a practical and user-friendly R-shiny tool to construct correlation networks from count tables. The details of data preprocessing, network construction, and visualization are automated, require no programming ability for the web version, and are highly customizable. Furthermore, SpeSpeNet detects clusters of co-occurring microbes and visualizes the relationship between these clusters or individual taxa and environmental metadata.

Find the full paper presenting SpeSpeNet here:

## Description files on this github:

app.R: Script for the ui and server of SpeSpeNet.

nw_functions.R: Script containing functions sourced by app.R.

Launch_tool.R: Script to load in libraries and launch the desktop version of the tool.

Manual.docx: In depth description of how to use SpeSpeNet.

Example_data: Map containing three data sets that are in SpeSpeNet compatible format and can be used to try SpeSpeNet. The three datasets correspond to the three case studies presented in the paper (Hauptfeld et al, 2022; Meisner et al, 2018; Brenzinger et al, 2021).

Figure_1.R: Rscript with linearized SpeSpeNet code used to produce Figure 1 of the paper.

Figure_2.R: Rscript with linearized SpeSpeNet code used to produce Figure 2 of the paper.

Figure3_FigureS1_FigureS2_FigureS4.R: Figure_1.R: Rscript with linearized SpeSpeNet code used to produce Figure 3, Figure S1, Figure S2 and Figure S4 of the paper.

## Remarks

Currently it is very difficult to obtain a workable installation of R to run SpeSpeNet due to the MGnifyR package no longer being compatible with the other packages. Hence the webtool is recommended.

## Contact

For issues, requests or discussions please open an issue on the github or contact [a.l.vaneijnatten\@uu.nl](mailto:a.l.vaneijnatten@uu.nl){.email}.
