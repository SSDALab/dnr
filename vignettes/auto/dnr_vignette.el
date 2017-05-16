(TeX-add-style-hook
 "dnr_vignette"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "12pt")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "amsmath"
    "amsthm"
    "amssymb"
    "geometry"
    "graphicx"
    "setspace"
    "verbatim")
   (TeX-add-symbols
    "myeq")
   (LaTeX-add-labels
    "sec:analysis-beach-data"
    "sec:model-fitting-beach"
    "sec:pred-beach-data"
    "sec:model-fixed-vertex"
    "sec:time-seri-param"
    "sec:performance-metrics"))
 :latex)

