# 3D Morphometry Toolbox

A Python-based ArcGIS toolbox to generate 3D morphometric data of elongated landscape features at user-defined transect-segmented intervals.

Journal: Computers and Geosciences

Author Details: Gwyneth Rivers*, Robert Storrar, Andrew Jones
* G.Rivers@shu.ac.uk; Sheffield Hallam University, Howard Street, Sheffield, S1 1WB


Programme Description:

The 3D Morphometry Toolbox is designed to efficiently extract and calculate 3D morphometric data of elongated landscape geomorphology. Specifically, the toolbox calculates detailed 3D morphometrics at user-defined transect-segmented intervals along the profile of a given landscape feature.

The toolbox comprises two scripts, each to be run as independent ArcGIS tools: a primary tool '3D Morphometry Tool' and a secondary tool 'Average Feature Morphometry'. The primary tool is calculates transect morphometry, whilst the secondary tool can be executed to average the transect morphometrics per parent feature, if desired.

The tools are written in Python [v2.7] and incorporate Python libraries from ‘ArcPy’ [ArcGIS 10.0-10.6], os, ‘Pandas’ [McKinney et al., 2010], and embeds Python code for tools, ‘Transect2.0’ [created by Mateus Vidotti Ferreira], and ‘Create Points on Lines’ [created by Ian Broad]. The toolbox needs to be downloaded and imported into the general ArcGIS ‘ArcToolbox’ workspace prior to use. 


Software Requirements and Availability:

The toolbox is intended to work within ArcGIS 10.1 [ArcMap; ESRI, 2018] and subsequent versions (including ArcGIS Pro [ESRI, 2020]). A ‘3D Analyst’ and ‘Spatial Analyst’ license is required. 

To use the toolbox, download and save the 3D Morphometry Toolbox zip file, open ArcGIS, import the toolbox into the general ArcGIS ArcToolbox and select the tool. Refer to the instruction guide and demonstration video for more detailed guidance if required.


Embedded Codes: 

Ian Broad. [Create Points on Lines Tool]. Web: www.ianbroad.com 

Mateus Vidotti Ferreira. [Transect2.0 Tool]. Email: mateusvidotti@yahoo.com.br 

Demonstration Video: 

https://youtu.be/9P4WkQxHbiU
