# LattUtil
Utility for handling VASP structure file

This is a utility program for handling VASP structure file (i.e. POSCAR).                                   
                                                                                                                
List of capabilities:                                                                                           
1. Find and display the symmetry of the system (info of space group, etc.)                                     
2. Regulate the structure with a fine precision                                                                
   (this is sometimes useful since the output of some software is only in single float prec).                  
3. Re-define the lattice (you may use this function to build a supercell)                                      
                                                                                                                
Note 1: LattUtil relies on spglib (later than 1.11) and armadillo,                                              
        you need to compile these two libs before compiling LattUtil.                                           
Note 2: LattUtil employs regular expression to handle text files,                                               
        therefore you need a relatively later version of gcc (later than 4.9).                                  
                                                                                                                
Developed by WSQ in May 2022.                                                                                   
