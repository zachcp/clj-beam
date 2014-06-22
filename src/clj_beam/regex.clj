(ns clj-beam.reagex
  (:require
   [clojure.set :as cset]
   [clojure.string :as string]
   [schema.core :as s]
   [clj-beam.schemas :as schemas]))


; example Strings to pick out
; OH3+ 2H 235U C@@H C@H I- Na+ Br- H+ Fe+2 OH- Fe++ OH3+ NH4+

; java regex reference here:
; http://www.tutorialspoint.com/java/java_regular_expressions.htm

(def examples ["OH3+" "Fe+++" "35Cl:classinfo" "27Fe+2" "2H" "235U" "C@@H" "C@H" "I-" "I-2" "Na+" "Br-" "H+" "Fe+2" "OH-" "Fe++" "OH3+" "NH4+"])

(defn bracketatominfo [a]
  "a bracketed   atom can may have the followinwg:
     isotope:    a nubmer before the symbol
     symbol:     the atomic symbol
     chirality:  for chiral centers
     protons:    H with an optional digit after
     charge:     multiple +/- with optional digit
     class:      optional field to specify extra information"
  (let [;get basic charge, proton, and isotope information
        element   (let  [ el (re-seq #"(?i)[abcdefgijklmnopqrstuvmwyz]+" a) ]
                    (if (nil? el)  "H" el)) ;pattern excludes "H" for ease. add it back if necessary
        isotope   (re-find #"^[0-9]+" a)
        chiral    (re-seq #"[@]+" a)
        chiralH   (re-seq #"[@]+[H]" a)
        charge    (first (re-seq #"[+-]+" a))
        chargenum (first (map read-string (re-find #"[+-](\d)" a)))
        protons   (re-seq #"[H]+" a)
        protnum   (second ( map read-string (re-find #"[H](\d)" a)))
        atomclass (re-seq #":(?i)(\w+)" a)

        ;resolve ambiguities
        fin_charge (cond
                     chargenum chargenum
                     charge  ( if (= (first charge) \+)
                                  (count charge)
                                  (- (count charge)))
                    :else 0)

        fin_protons (cond
                     protnum protnum
                     protons 1
                    :else 0)
        ]

    ;[a element isotope chiral chiralH charge chargenum fin_charge protons protnum atomclass]


     { :element   (first element)
       :isotope   isotope
       :charge    fin_charge
       :hydrogens fin_protons
       :atomclass atomclass }

;    [a protons protnum]

    ))


(map bracketatominfo examples)