# Thesis
## Compilation
Compile with _pdfLaTeX_ and _biber_ running this sequence of commands:
```
$ pdflatex main.tex
$ biber main
$ pdflatex main.tex
$ pdflatex main.tex
```
List of required packages is in preamb.tex.

For an automated way to install packages, bibliograhy generation and one line compilation, try using _tectonic_:
```
$ tectonic main.tex
```
You will still need to install biber manually, but manual tex-live installation is not required.
Instruction on how to get it on its [github page](https://github.com/tectonic-typesetting/tectonic).