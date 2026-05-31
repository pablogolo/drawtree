DRAWTREE - A small tree-drawing program.

First argument is name of file to read; if parenthetical tree uses only numbers (not taxon names), you can list the taxon names in a separate file, given to DRAWTREE as second argument.

Can compile with 

    gcc -O3 -o drawtree drawtree.c

In Windows, compiling in DOS mode produces a nice-looking tree. You can use it to quickly inspect TNT tree files. 

In other OSs, you may have to set the screen font to Tred (comes with TNT, included here as well). Otherwise, you can compile the code to draw the "poorman" tree, with

    gcc -DPOORMAN -o drawtree drawtree.c
    
