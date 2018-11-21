
DEW_2_ver converts thermodynamic data from the Harrison & Sverjensky DEW Excel spreadsheet (www.dewcommunity.org/resources.html) 
to Perple_X format/units. 

DEW_2_ver reads two input files (examples should accompany this README file):

   DEW_input.txt - the file containing the DEW format data.
   DEW_elements_input.dat - a file defining the elements to be retained and conversion options.

DEW_2_ver writes a file named DEW_2_ver_output.dat containing the Perple_X formatted data, this file does not
have the standard Perple_X thermodynamic data file header. Therefore the header must be copied from another 
thermodynamic data file, or the data in DEW_2_ver_output.dat must be added to an existing thermodynamic data 
file (e.g., one containing condensed and solvent phase data).

1) To prepare the DEW_input.txt:

    a) Open the DEW spreadsheet in Excel
    b) If desired, change the data options on the "Aqueous Species Options" page.
    c) Change to the "Aqueous Species Table", verify the columns and units of the DEW data correspond to:

       name | formula |   Gf      |  Hf       | So            | Vo        | Cpo           | a1 x 10         | a2 x 10-2  | a3                | a4 x 10-4   | c1            | c2 x 10-4   | omega x 10-5  | Z | Comments
                      | cal mol-1 | cal mol-1 | cal mol-1 K-1 | cm3 mol-1 | cal mol-1 K-1 | cal mol-1 bar-1 | cal mol-1  | cal K mol-1 bar-1 | cal K mol-1 | cal mol-1 K-1 | cal K mol-1 | cal mol-1     |   |

    d) Copy the data, without the header and the blank first column, and paste it into a text file named
       DEW_input.txt, replace any tab characters with blanks and save the file (alternatively use the "Save 
       as" option to save the data to a space delimited text file, in this case it is probably necessary to
       change the file-type suffix from ".prn" to ".txt", also it is wise to verify that species names are
       left justified).
    e) BEWARE!!! a few entries in the DEW spreadsheet contain no enthalpy value, DEW_2_VER recognizes this if it
       reads an alphabetic character for an entry after reading only 12 of the 13 numeric values normally
       present for each species. This method will fail if no comment is present or the comment begins with 
       a numeric character (e.g., as in a date). To avoid this problem, scan the "Aqueous Species Table" and 
       either add a number for any missing enthalpy value or insert a non-numeric character in front of any
       comment that begins with a number. 

2) To modify/prepare DEW_elements_input.dat see the comments in the example included herewith. 

3) Both input files and the DEW_2_ver program must be located in the same directory/folder.

4) Run DEW_2_ver. The program writes three types of diagnostic to the user console:

   a) "truncated name from HEXANE,AQ to HEXANE,A" - indicates that the DEW species name has been truncated
      to the 8 character format required by Perple_X. This diagnostic may be important if truncation 
      results in two different species having the same 8 name in Perple_X. Perple_X will identify such
      cases as replicates and terminate with an error message. To avoid this modify the names in the 
      Perple_X formatted output. 
   b) "probably missing enthalpy for:KCl,aq" - indicates the DEW entry for KCl,aq has only 12 entries, 
      invariably this has been the enthalpy, which is not required by DEW_2_ver unless the HSC 
      convention flag is set to .true. in DEW_elements_input.dat. 
   c)  "found a repeat group these seem to be useless will ignore, CHECK THIS FORMULA!
        mugga: Na(CH3COO)2(-1)" - formulae written in the DEW spreadsheet rarely employ repeat groups
       (e.g., ...(OH)2) and when they are present they are potentially ambiguous with the (superfluous)
       charge notation. For these reasons, DEW_2_ver does not attempt to interpret possible repeat 
       group stoichiometry. If the number following a parenthesis in the DEW formula is really a 
       stoichiometric coefficient, then the formula in DEW_2_ver_output.dat must be corrected. In the
       current (May 19, 2017) version of the data this is true for 3 species: B(OH)3(0), Na(CH3COO)2(-1),
       and Fe(CH3COO)2(0) [yes, I know, the stoichiometric coeffiecient is not in brackets, but the 
       charge is. I am not going to program for only 3 instances.]

