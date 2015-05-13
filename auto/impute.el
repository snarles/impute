(TeX-add-style-hook "impute"
 (lambda ()
    (LaTeX-add-bibliographies
     "gwas"
     "books")
    (LaTeX-add-labels
     "domcond")
    (TeX-add-symbols
     "tr"
     "E"
     "diag"
     "argmax"
     "Cov")
    (TeX-run-style-hooks
     "mathtools"
     "graphicx"
     "fancyhdr"
     "color"
     "ifthen"
     "array"
     "epsfig"
     "amsmath"
     "amssymb"
     "latex2e"
     "art12"
     "article"
     "12pt")))

