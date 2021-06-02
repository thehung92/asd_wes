#!/usr/bin/env Rscript
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
install.packages('pacman', dependencies=TRUE)
pacman::p_install("CMplot")
pacman::p_install("Gviz")
pacman::p_install("lemon")
pacman::p_install("biomaRt")
pacman::p_install("grid")
pacman::p_install("gridtext")
pacman::p_install("ggpubr")
pacman::p_install("gprofiler2")
pacman::p_install("openxlsx")