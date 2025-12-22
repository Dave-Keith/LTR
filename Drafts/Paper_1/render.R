#Step 1: render the csasdown project
# (working directory must be set to Framework/DataInputs/csasdown)
setwd("D:/Github/ICM/Drafts/Paper_1/")
csasdown::render()
# Step 2: Add initial pages of a typical resdoc (title page, etc.)
resdoc <- "_book/paper_1.docx"

print(x, target = "_book/resdoc-english.docx")

# Add appendix with diff numbering
rmarkdown::render("appendix.Rmd")
supplement <- "appendix.docx"
resdoc <- officer::read_docx(path = "_book/paper_1.docx")
x <- officer::body_add_docx(resdoc, supplement)
# Step 4: Save
# If error message, look for temp file in RStudio Files panel. Select box, and delete in Rstudio. Temp file may not be visible in File explorer.
print(x, target = "_book/final_draft.docx")
