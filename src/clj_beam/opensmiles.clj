; Reference for Implementing a Parser
; From OpenSmiles.org
; http://www.opensmiles.org/opensmiles.html


;Atoms
;atom ::= bracket_atom aliphatic_organic aromatic_organic "*"


;Organic Subset
(def organicsubset {:aliphatic_organic ["B" "C" "N" "O" "S" "P" "F" "Cl" "Br" "I"]
                    :aromatic_organic  ["b" "c" "n" "o" "s" "p"]})

;Bracket Atoms
(def bracket_atom  [ "\[" :isotope? :symbol :chiral? :hcount? :charge? :class? "\]" ] )

(def  symbol1 [ element_symbols aromatic_symbols "*"] )
;isotope ::= NUMBER

(def element_symbols ["H" "He" "Li" "Be" "B" "C" "N" "O" "F" "Ne" "Na" "Mg" "Al" "Si"
                      "P" "S" "Cl" "Ar" "K" "Ca" "Sc" "Ti" "V" "Cr" "Mn" "Fe" "Co" "Ni"
                      "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr" "Rb" "Sr" "Y" "Zr" "Nb" "Mo" "
                      Tc" "Ru" "Rh" "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "I" "Xe" "Cs" "Ba" "Hf"
                      "Ta" "W" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl" "Pb" "Bi" "Po" "At" "Rn" "Fr"
                      "Ra" "Rf" "Db" "Sg" "Bh" "Hs" "Mt" "Ds" "Rg" "Cn" "Fl" "Lv" "La" "Ce"
                      "Pr" "Nd" "Pm" "Sm" "Eu" "Gd" "Tb" "Dy" "Ho" "Er" "Tm" "Yb" "Lu" "Ac" "Th"
                      "Pa" "U" "Np" "Pu" "Am" "Cm" "Bk" "Cf" "Es" "Fm" "Md" "No" "Lr"])

(def aromatic_symbols ["b" "c" "n" "o" "p" "s" "se" "as" ])


;Chirality
(def chiral ["@" "@@"
             "@TH1" "@TH2"
             "@AL1" "@AL2"
             "@SP1" "@SP2" "@SP3"
             "@TB1" "@TB2" "@TB3" "@TB4" "@TB5" "@TB6" "@TB7" "@TB8" "@TB9" "@TB10" "@TB11" "@TB12" "@TB13" "@TB14" "@TB15" "@TB16" "@TB17" "@TB18" "@TB19" "@TB20"
             "@OH1" "@OH2" "@OH3" "@OH4" "@OH5" "@OH6" "@OH7" "@OH8" "@OH9" "@OH10" "@OH11" "@OH12" "@OH13" "@OH14" "@OH15" "@OH16" "@OH17" "@OH18" "@OH19" "@OH20" "@OH21" "@OH22" "@OH23" "@OH24" "@OH25" "@OH26" "@OH27" "@OH28" "@OH29" "@OH30"
             ;"@TB" DIGIT DIGIT
             ;"@OH" DIGIT DIGIT
             ]


;Hydrogens
;hcount ::= "H" "H" DIGIT

;Charge
;charge ::= "-" "-" DIGIT? DIGIT "+" "+" DIGIT? DIGIT "--" deprecated "++" deprecated

;Atom Class
;class ::= ":" NUMBER

;Bonds

(def bond [ "-" "=" "#" "$" ":" "\/" "\\" ])
;ringbond ::= bond? DIGIT bond? "%" DIGIT DIGIT
;branched_atom ::= atom ringbond* branch*
;branch ::= "(" chain ")" "(" bond chain ")" "(" dot chain ")"
;chain ::= branched_atom chain branched_atom chain bond branched_atom chain dot branched_atom
;dot ::= "."

;SMILES STRINGS
;smiles ::= terminator chain terminator
;terminator ::= SPACE TAB LINEFEED CARRIAGE_RETURN END_OF_STRING