# scStudio
A User-Friendly Web Application Empowering Non-Computational Users with Intuitive scRNA-seq Data Analysis

scStudio is a user-friendly, comprehensive and modular web-based application designed to democratize scRNA-seq data analysis. scStudio is equipped with a suite of features designed to streamline data retrieval and analysis with both flexibility and ease, including automated dataset retrieval from the Gene Expression Omnibus (GEO). Users can also upload their own datasets in a variety of formats, integrate multiple datasets, and tailor their analyses using a wide range of flexible methods with options for parameter optimization. The application supports all the essential steps required for scRNA-seq data analysis, including in-depth quality control, normalization, dimensionality reduction, clustering, differential expression and functional enrichment analysis. scStudio also tracks the history of analyses, supports session data storage and export, and facilitates collaboration through data sharing features. By developing scStudio as a user-friendly interface and scalable architecture, we address the evolving needs of scRNA-seq research, making advanced data analysis accessible and manageable while accommodating future developments in the field.

# Development

scStudio was developed as a web-based application using the R Shiny package (version 1.10.0). Each main analysis step is implemented as a separate Shiny app. This modular approach allows users to run multiple types of analyses simultaneously, each in its own dedicated app window. As different methods complete their computations, the results are automatically updated, ensuring that users have real-time access to their data and analyses. 

To begin using scStudio, start by creating a working session in the Upload Data interface. This interface allows you to:

*Retrieve datasets directly from GEO

*Upload your own data in various formats

*Load public working sessions from previously processed datasets

Each session is assigned a unique token, which you can use to save and resume your analysis, share your session with others, or download your data files for local use.
[IMPORTANT:] Sessions that remain inactive for over two weeks will be automatically removed from the server. If you wish to preserve your session for a longer period, please contact: ana.bica@gimm.pt.

# Upload data interface

## Gene Expression Omnibus interface

scStudio enables the automatic retrieval of scRNA-seq data from GEO. To get started, simply enter the GEO accession ID of your dataset and click “Load.” Once the data is loaded, you can select the samples of interest and click “Add samples” to generate a unified count matrix for downstream analysis.

Note: Currently supported data formats include tabular files (CSV, TSV, TXT), Excel spreadsheets (XLSX), Cell Ranger Market Exchange Format (MEX) files commonly used for 10X Genomics data, and HDF5 (.h5) files.
   













