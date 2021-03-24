rmarkdown::render("Scripts/01_maln.Rmd", output_file = "01_maln_pre2000.md", params = list(restrict = "pre"))
rmarkdown::render("Scripts/01_maln.Rmd", output_file = "01_maln_pos2000.md", params = list(restrict = "post"))
rmarkdown::render("Scripts/01_maln.Rmd", output_file = "01_maln.md", params = list(restrict = "all"))
