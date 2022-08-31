# Thesis
## Compilation
Recommended way to compile is to use _tectonic_, which gives an automated way to install packages, bibliograhy generation (using biber) and one line compilation:
```
$ tectonic main.tex
```
Instructions on how to get it from its [github page](https://github.com/tectonic-typesetting/tectonic).
Manual installation of biber is required:
```
$ sudo apt install biber
```

Alternatively, you should be able to compile with _pdfLaTeX_ and _biber_ running this sequence of commands:
```
$ pdflatex main.tex
$ biber main
$ pdflatex main.tex
$ pdflatex main.tex
```
List of required packages is in preamb.tex.